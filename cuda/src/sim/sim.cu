#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdbool.h>

#include <io/error.h>
#include <io/error_cu.h>
#include <sim/kernel.h>
#include <sim/sim.h>

static inline FluidVoxel_t* index_voxel(SimState_t* sim, int x, int y)
{
    return &(sim->voxels[y * (sim->params.dims.x) + x]);
}

FluidsimError_t _simStateAllocOnDevice(SimState_t** d_onDevice, SimParams_t params){
    size_t arrSize = params.dims.x * params.dims.y * sizeof(FluidVoxel_t);
    FluidVoxel_t *arr, *arr2;
    cudaMalloc(&arr, arrSize);
    cudaMalloc(&arr2, arrSize);

    SimState_t h_init = {
        .frame = 0,
        .params = params,
        .voxels = arr,
        .d_deviceStatePtr = NULL,
        ._d_voxels_old = arr2
    };

    cudaMalloc(d_onDevice, sizeof(SimState_t));
    cudaMemcpy(*d_onDevice, &h_init, sizeof(SimState_t), cudaMemcpyHostToDevice);

    return FSE_OK;
}

FluidsimError_t _simStateCopyToDevice(SimState_t* d_onDevice, SimState_t* h_onHost){
    // Get array pointers
    // N.B. - h_onDevice will be a struct ON THE HOST containing information on
    //        the simulation, as well as a pointer to the array ON THE DEVICE
    //        where voxels are stored.
    SimState_t h_onDevice;
    cudaMemcpy(&h_onDevice, d_onDevice, sizeof(SimState_t), cudaMemcpyDeviceToHost);

    // Copy fluid voxels
    if(
        h_onDevice.params.dims.x * h_onDevice.params.dims.y !=
         h_onHost->params.dims.x *  h_onHost->params.dims.y
    ){
        // Simulation size change
        fprintf(
            stderr,
            "ERROR: Simulation size change from (%d x %d) to (%d x %d). "
            "Simulation state must be re-initialized if the size changes!\n",
            h_onHost->params.dims.x, h_onHost->params.dims.y, 
            h_onDevice.params.dims.x, h_onDevice.params.dims.y
        );
        return FSE_BAD_SIM_STATE;
    }
    size_t arrayLen = h_onHost->params.dims.x * 
                      h_onHost->params.dims.y * 
                      sizeof(FluidVoxel_t);
    cudaMemcpy(h_onDevice.voxels, h_onHost->voxels, arrayLen, cudaMemcpyHostToDevice);
    cudaMemcpy(h_onDevice._d_voxels_old, h_onHost->voxels, arrayLen, cudaMemcpyHostToDevice);

    // Copy base struct, substituting our new voxel array and leaving a few
    // other values untouched
    SimState_t h_toCopy = {
        .frame = h_onHost->frame,
        .params = h_onHost->params,
        .voxels = h_onDevice.voxels,
        .d_deviceStatePtr = NULL,
        ._d_voxels_old = h_onDevice._d_voxels_old
    };
    cudaMemcpy(d_onDevice, &h_toCopy, sizeof(SimState_t), cudaMemcpyHostToDevice);

    return FSE_OK;
}

FluidsimError_t _simStateCopyToHost(SimState_t* d_onDevice, SimState_t* h_onHost){
    // Preserve host voxel-array address and device-state pointer
    FluidVoxel_t* h_arr = h_onHost->voxels;
    SimState_t* d_hostDevPtr = h_onHost->d_deviceStatePtr;

    // Copy base struct
    cudaMemcpy(h_onHost, d_onDevice, sizeof(SimState_t), cudaMemcpyDeviceToHost);

    // Preserve device voxel-array
    FluidVoxel_t* d_arr = h_onHost->voxels;

    // Restore host voxel-array address and device-state pointer
    h_onHost->voxels = h_arr;
    h_onHost->d_deviceStatePtr = d_hostDevPtr;

    // Copy voxel array
    size_t arrayLen = h_onHost->params.dims.x * 
                      h_onHost->params.dims.y * 
                      sizeof(FluidVoxel_t);
    cudaMemcpy(h_arr, d_arr, arrayLen, cudaMemcpyDeviceToHost);

    return FSE_OK;
}

FluidsimError_t _simStateFreeOnDevice(SimState_t* d_onDevice){
    // Retrieve data (we need the array pointer from here)
    SimState_t h_onDevice;
    cudaMemcpy(&h_onDevice, d_onDevice, sizeof(SimState_t), cudaMemcpyDeviceToHost);

    // Free voxel array
    cudaFree(h_onDevice.voxels);

    // Free the struct itself
    cudaFree(d_onDevice);
    return FSE_OK;
}

FluidsimError_t initSimState(SimParams_t params, SimState_t* h_toInit){
    // Initialize locally
    h_toInit->params = params;
    h_toInit->frame = 0;

    int arrSize = params.dims.x*params.dims.y*sizeof(FluidVoxel_t);
    
    FluidVoxel_t* voxels = (FluidVoxel_t*) malloc(arrSize);

    if(!voxels){
        return FSE_MALLOC_FAILED;
    }

    h_toInit->voxels = voxels;

    // Initialize on device
    SimState_t* d_onDevice = NULL;
    FluidsimError_t err = _simStateAllocOnDevice(&d_onDevice, params);
    if(err != FSE_OK) return err;

    // Copy to device
    err = _simStateCopyToDevice(d_onDevice, h_toInit);
    if(err != FSE_OK) return err;

    // Set device pointer
    h_toInit->d_deviceStatePtr = d_onDevice;

    // Set up CUDA environment
    dim3 dimGrid = {
        (unsigned int)((h_toInit->params.dims.x + TILE_WIDTH - 1) / TILE_WIDTH), 
        (unsigned int)((h_toInit->params.dims.y + TILE_WIDTH - 1) / TILE_WIDTH)
    };
    dim3 dimBlock = {
        TILE_WIDTH,
        TILE_WIDTH
    };

    // Launch initialization kernel
    InitializerKernel<<<dimGrid, dimBlock>>>(d_onDevice);

    // Check for errors
    fseCuChk(cudaPeekAtLastError());
    fseCuChk(cudaDeviceSynchronize());

    // Copy-back to host
    err = _simStateCopyToHost(d_onDevice, h_toInit);

    return err;
}

FluidsimError_t destroySimState(SimState_t* toDestroy){
    // Avoid duplicate frees
    if(toDestroy->d_deviceStatePtr == NULL) return FSE_NULL_PTR;

    // Free device memory (the entire struct)
    FluidsimError_t res = _simStateFreeOnDevice(toDestroy->d_deviceStatePtr);
    toDestroy->d_deviceStatePtr = NULL;

    // Free local memory (just the voxel array)
    free(toDestroy->voxels);

    if(res != FSE_OK) return res;

    return FSE_OK;
}

FluidsimError_t syncSimStateToHost(SimState_t* h_onHost){
    return _simStateCopyToHost(h_onHost->d_deviceStatePtr, h_onHost);
}

FluidsimError_t syncSimStateToDevice(SimState_t* h_onHost){
    return _simStateCopyToDevice(h_onHost->d_deviceStatePtr, h_onHost);
}

FluidsimError_t doFrame(KernelSet_t kernel, SimState_t* sim, float* time){
    // Set up CUDA environment
    int trueTileWidth = kernel.tileWidth - 2*kernel.tileOverlap;
    int trueTileHeight = kernel.tileHeight - 2*kernel.tileOverlap;
    dim3 dimGrid = {
        (unsigned int)((sim->params.dims.x + trueTileWidth - 1) / trueTileWidth), 
        (unsigned int)((sim->params.dims.y + trueTileHeight - 1) / trueTileHeight)
    };
    dim3 dimBlock = {
        kernel.tileWidth,
        kernel.tileHeight
    };
    cudaEvent_t start, stop;
    fseCuChk(cudaEventCreate(&start));
    fseCuChk(cudaEventCreate(&stop));

    // Launch kernel
    // kernel<<<dimGrid, dimBlock>>>(sim->d_deviceStatePtr);

    if(kernel.full){
        // Run
        fseCuChk(cudaEventRecord(start));
        kernel.full<<<dimGrid, dimBlock>>>(sim->d_deviceStatePtr);
        fseCuChk(cudaEventRecord(stop));
        // Check for errors
        fseCuChk(cudaPeekAtLastError());
        fseCuChk(cudaDeviceSynchronize());

    }else{
        // Collide
        fseCuChk(cudaEventRecord(start));
        kernel.collide<<<dimGrid, dimBlock>>>(sim->d_deviceStatePtr);
        // Check for errors
        fseCuChk(cudaPeekAtLastError());
        fseCuChk(cudaDeviceSynchronize());

        // Exchange
        kernel.exchange<<<dimGrid, dimBlock>>>(sim->d_deviceStatePtr);
        // Check for errors
        fseCuChk(cudaPeekAtLastError());
        fseCuChk(cudaDeviceSynchronize());

        // Stream
        kernel.stream<<<dimGrid, dimBlock>>>(sim->d_deviceStatePtr);
        // Check for errors
        fseCuChk(cudaPeekAtLastError());
        fseCuChk(cudaDeviceSynchronize());

        // Bounceback
        kernel.bounceBack<<<dimGrid, dimBlock>>>(sim->d_deviceStatePtr);
        fseCuChk(cudaEventRecord(stop));
        // Check for errors
        fseCuChk(cudaPeekAtLastError());
        fseCuChk(cudaDeviceSynchronize());
    }

    fseCuChk(cudaEventSynchronize(stop));

    if(time != NULL){
        fseCuChk(cudaEventElapsedTime(time, start, stop));
    }


    return FSE_OK;
}

FluidsimError_t getCurl(SimState_t* sim, float* curl){
    if((!curl) || (!sim)){
        return FSE_NULL_PTR;
    }

    for(int y = 1; y < sim->params.dims.y - 1; ++y){
        for(int x = 1; x < sim->params.dims.x - 1; ++x){
            int dex = y*sim->params.dims.y + x;
            curl[dex] = index_voxel(sim, x + 1, y)->velocity.y
                      + index_voxel(sim, x - 1, y)->velocity.y
                      + index_voxel(sim, x, y + 1)->velocity.x
                      + index_voxel(sim, x, y - 1)->velocity.x;
        }
    }

    return FSE_OK;
}

FluidsimError_t initLogFile(FILE* f, SimParams_t p){
    if(!f){
        return FSE_NULL_PTR;
    }

    fprintf(
        f,
        "curl: f[%5i, %5i];\r\n"
        "density: f[%5i, %5i];\r\n"
        "speed: f[%5i, %5i];\r\n"
        "barriers: ?[%5i, %5i];\r\n"
        "FSIM",
        p.dims.x, p.dims.y,
        p.dims.x, p.dims.y,
        p.dims.x, p.dims.y,
        p.dims.x, p.dims.y
    );

    return FSE_OK;
}

FluidsimError_t writeLogFrame(FILE* f, SimState_t* sim){
    size_t numEntries = sim->params.dims.x * sim->params.dims.y;
    float* curls = (float*)malloc(numEntries * sizeof(float));
    float* densities = (float*)malloc(numEntries * sizeof(float));
    float* speeds = (float*)malloc(numEntries * sizeof(float));
    uint8_t* barriers = (uint8_t*)malloc(numEntries * sizeof(uint8_t));

    FluidsimError_t e = getCurl(sim, curls);
    if(e != FSE_OK) return e;

    int dex = 0;
    for(int y = 0; y < sim->params.dims.y; ++y){
        for(int x = 0; x < sim->params.dims.x; ++x){
            FluidVoxel_t* v = index_voxel(sim, x, y);
            densities[dex] = v->density;
            speeds[dex] = sqrt((v->velocity.x * v->velocity.x) + (v->velocity.y * v->velocity.y));
            barriers[dex] = v->is_barrier;
            ++dex;
            // printf("[ (%d,%d) v=<%f,%f> rho=%f barrier=%s ]\n", x, y, v->velocity.x, v->velocity.y, v->density, v->is_barrier ? "yes" : "no");
        }
    }


    if(!fwrite(curls,     sizeof(float),   numEntries, f)) return FSE_FILE_IO_FAILURE;
    if(!fwrite(densities, sizeof(float),   numEntries, f)) return FSE_FILE_IO_FAILURE;
    if(!fwrite(speeds,    sizeof(float),   numEntries, f)) return FSE_FILE_IO_FAILURE;
    if(!fwrite(barriers,  sizeof(uint8_t), numEntries, f)) return FSE_FILE_IO_FAILURE;

    return FSE_OK;
}
