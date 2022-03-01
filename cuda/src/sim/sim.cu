#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdbool.h>

#include <io/error.h>
#include <sim/kernel.h>
#include <sim/sim.h>

static inline FluidVoxel_t* index_voxel(SimState_t* sim, int x, int y)
{
    return &(sim->voxels[y * (sim->params.dims.x) + x]);
}

FluidsimError_t _simStateAllocOnDevice(SimState_t** d_onDevice, SimParams_t params){
    size_t arrSize = params.dims.x * params.dims.y * sizeof(FluidVoxel_t);
    FluidVoxel_t* arr;
    cudaMalloc(&arr, arrSize);

    SimState_t h_init = {
        .frame = 0,
        .params = params,
        .voxels = arr
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
            "Simulation state must be re-allocated if the size changes!\n",
            h_onHost->params.dims.x, h_onHost->params.dims.y, 
            h_onDevice.params.dims.x, h_onDevice.params.dims.y
        );
        return FSE_BAD_SIM_STATE;
    }
    size_t arrayLen = h_onHost->params.dims.x * 
                      h_onHost->params.dims.y * 
                      sizeof(FluidVoxel_t);
    cudaMemcpy(h_onDevice.voxels, h_onHost->voxels, arrayLen, cudaMemcpyHostToDevice);

    // Copy base struct, substituting our new voxel array
    SimState_t h_toCopy = {
        .frame = h_onHost->frame,
        .params = h_onHost->params,
        .voxels = h_onDevice.voxels
    };
    cudaMemcpy(d_onDevice, &h_toCopy, sizeof(SimState_t), cudaMemcpyHostToDevice);

    return FSE_OK;
}

FluidsimError_t _simStateCopyToHost(SimState_t* d_onDevice, SimState_t* h_onHost){
    // Preserve host voxel-array address
    FluidVoxel_t* h_arr = h_onHost->voxels;

    // Copy base struct
    cudaMemcpy(h_onHost, d_onDevice, sizeof(SimState_t), cudaMemcpyDeviceToHost);

    // Preserve device voxel-array address
    FluidVoxel_t* d_arr = h_onHost->voxels;

    // Restore host voxel-array address
    h_onHost->voxels = h_arr;

    // Copy voxel array
    size_t arrayLen = h_onHost->params.dims.x * 
                      h_onHost->params.dims.y * 
                      sizeof(FluidVoxel_t);
    cudaMemcpy(h_arr, d_arr, arrayLen, cudaMemcpyDeviceToHost);

    return FSE_OK;
}

FluidsimError_t initSimState(SimParams_t params, SimState_t* toInit){
    toInit->params = params;
    toInit->frame = 0;

    int arrSize = params.dims.x*params.dims.y*sizeof(FluidVoxel_t);
    
    FluidVoxel_t* voxels = (FluidVoxel_t*) malloc(arrSize);

    if(!voxels){
        return FSE_MALLOC_FAILED;
    }

    toInit->voxels = voxels;
    return FSE_OK;
}

FluidsimError_t destroySimState(SimState_t* toDestroy){
    free(toDestroy->voxels);
    return FSE_OK;
}

FluidsimError_t doFrame(Kernel_t kernel, SimState_t* sim){
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
        }
    }

    if(!fwrite(curls,     sizeof(float),   numEntries, f)) return FSE_FILE_IO_FAILURE;
    if(!fwrite(densities, sizeof(float),   numEntries, f)) return FSE_FILE_IO_FAILURE;
    if(!fwrite(speeds,    sizeof(float),   numEntries, f)) return FSE_FILE_IO_FAILURE;
    if(!fwrite(barriers,  sizeof(uint8_t), numEntries, f)) return FSE_FILE_IO_FAILURE;

    return FSE_OK;
}
