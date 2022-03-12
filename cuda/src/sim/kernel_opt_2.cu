#include <io/error.h>
#include <sim/sim_struct.h>
#include <sim/kernel.h>

#include <sim/_kernel_util.h>

#include <stdio.h>
#include <stdbool.h>

#undef TILE_WIDTH
#define TILE_WIDTH 32
#define TILE_HEIGHT 8
#define TILE_OVERLAP 1

#define IM_X ((signed int)blockIdx.x * ((signed int)blockDim.x-2*TILE_OVERLAP) + (signed int)threadIdx.x - TILE_OVERLAP)
#define IM_Y ((signed int)blockIdx.y * ((signed int)blockDim.y-2*TILE_OVERLAP) + (signed int)threadIdx.y - TILE_OVERLAP)
#define IM_BASE_X ((signed int)blockIdx.x * ((signed int)blockDim.x-2*TILE_OVERLAP) - TILE_OVERLAP)
#define IM_OFFS(sim)(_VOXEL_IDX(sim, IM_X, IM_Y))

#define MY_TILE_OFFS (threadIdx.y * blockDim.x + threadIdx.x)

static inline __device__ FluidVoxel_t* _index_tile_voxel(FluidVoxel_t* tile, int xOffs, int yOffs){
    int offs = (threadIdx.y + yOffs)*blockDim.x + threadIdx.x + xOffs;
    return tile + offs;
}

static __device__ void _collide(SimState_t* state, int myX, int myY, FluidVoxel_t* oldVA){
    if(myX < state->params.dims.x && myY < state->params.dims.y && myX >= 0 && myY >= 0){
        const float omega = 1.0f / ((3.0f*state->params.viscosity) + 0.5f);

        FluidVoxel_t* myVoxel = _index_tile_voxel(oldVA, 0, 0);

        // if(!isShared(myVoxel)){
        //     printf("Tile coordinates (%d,%d), Pointer: %x\n")
        // }

        // printf("Pointer: %x\n", myVoxel);

        // float f = myVoxel->density;

        const static float ux_dirs[] = {0,-1,-1,-1, 0, 1, 1, 1, 0};
        const static float uy_dirs[] = {1, 1, 0,-1,-1,-1, 0, 1, 0};
        float rho = 0.0f, ux = 0.0f, uy = 0.0f;
        for(int i = 0; i < NUM_LATTICE_VECTORS; ++i){
            float lv = myVoxel->lattice_vectors.sequence[i];
            rho += lv;
            ux  += lv * ux_dirs[i];
            uy  += lv * uy_dirs[i];
        }

        ux /= rho;
        uy /= rho;

        myVoxel->density = rho;
        myVoxel->velocity.x = ux;
        myVoxel->velocity.y = uy;

        _mod_voxel_lattice_vectors(myVoxel, ux, uy, rho, omega);
    }
}

static __device__ void _stream(SimState_t* state, int myX, int myY, FluidVoxel_t* oldVA, FluidVoxel_t* newVA){
    if(
        myX < state->params.dims.x && myY < state->params.dims.y && myY >= 0 && myX >= 0 // Check image dimensions
        && threadIdx.x > 0 && threadIdx.x < blockDim.x - 1 // Exclude edge voxels (x)
        && threadIdx.y > 0 && threadIdx.y < blockDim.y - 1 // Exclude edge voxels (y)
    ){
        FluidVoxel_t* myVoxel = _index_tile_voxel(newVA, 0, 0);
        FluidVoxel_t* myOldVoxel = _index_tile_voxel(oldVA, 0, 0);

        myVoxel->density    = myOldVoxel->density;
        myVoxel->velocity.x = myOldVoxel->velocity.x;
        myVoxel->velocity.y = myOldVoxel->velocity.y;

        // ========= Copy lattice vectors =========
        // Copy from voxel to the south
        myVoxel->lattice_vectors.named.north 
            = _index_tile_voxel(oldVA, 0, -1)
                -> lattice_vectors.named.north;

        // --> South
        // Copy from the voxel to the north
        myVoxel->lattice_vectors.named.south
            = _index_tile_voxel(oldVA, 0, 1)
                -> lattice_vectors.named.south;

        // --> East
        // Copy from the voxel to the west
        myVoxel->lattice_vectors.named.east
            = _index_tile_voxel(oldVA, -1, 0)
                -> lattice_vectors.named.east;

        // --> West
        // Copy from the voxel to the east
        myVoxel->lattice_vectors.named.west
            = _index_tile_voxel(oldVA, 1, 0)
                -> lattice_vectors.named.west;

        // --> Northeast
        // Copy from the voxel to the southwest
        myVoxel->lattice_vectors.named.northeast
            = _index_tile_voxel(oldVA, -1, -1)
                -> lattice_vectors.named.northeast;

        // --> Northwest
        // Copy from the voxel to the southeast
        myVoxel->lattice_vectors.named.northwest
            = _index_tile_voxel(oldVA, 1, -1)
                -> lattice_vectors.named.northwest;

        // --> Southeast
        // Copy from the voxel to the northwest
        myVoxel->lattice_vectors.named.southeast
            = _index_tile_voxel(oldVA, -1, 1)
                -> lattice_vectors.named.southeast;

        // --> Southwest
        // Copy from the voxel to the northeast
        myVoxel->lattice_vectors.named.southwest
            = _index_tile_voxel(oldVA, 1, 1)
                -> lattice_vectors.named.southwest;
    }
}

static __device__ void _barrierBounceBack(SimState_t* state, int myX, int myY, FluidVoxel_t* newVA){
    if(        
        myX < state->params.dims.x && myY < state->params.dims.y && myY >= 0 && myX >= 0 // Check image dimensions
        && threadIdx.x > 0 && threadIdx.x < blockDim.x - 1 // Exclude edge voxels (x)
        && threadIdx.y > 0 && threadIdx.y < blockDim.y - 1 // Exclude edge voxels (y)
    ){
        FluidVoxel_t* myVoxel = _index_tile_voxel(newVA, 0, 0);
        if(myVoxel->is_barrier){
            // For every lattice vector except the zero vector
            for(int i = 0; i < LV_IM; ++i){
                // Find the voxel in the opposite direction (i.e. if i == LV_S,
                // get the voxel to the north)
                IntPoint_t voxelDelta = _voxel_delta_in_direction(LV_OPPOSITE_DIR_OF(i));

                // Bounce back
                _index_tile_voxel(newVA, voxelDelta.x, voxelDelta.y)
                    -> lattice_vectors.sequence[LV_OPPOSITE_DIR_OF(i)]
                        = myVoxel->lattice_vectors.sequence[i];
            }
        }
    }

    if(myX == 0 && myY == 0){
        ++(state->frame);
    }
}

static __device__ void _parl_copy(void* dst, void* src, size_t len, size_t numThreads, size_t myTid){
    // int passCount = (len + numThreads - 1) / numThreads;
    int passCount = len / numThreads;

    for(int i = 0; i < passCount; ++i){
        int myOffs = i*numThreads + myTid;
        char x = ((char*)src)[myOffs];
        ((char*)dst)[myOffs] = x;
    }

    int myOffs = passCount*numThreads + myTid;
    if(myOffs < len){
        char x = ((char*)src)[myOffs];
        ((char*)dst)[myOffs] = x;
    }
}

static __device__ void _populateRow(SimState_t* state, FluidVoxel_t* rowPtr){

    if(IM_Y >= 0 && IM_Y < state->params.dims.y){
        FluidVoxel_t* dstPtr = rowPtr;
        FluidVoxel_t* srcPtr = state->_d_voxels_old + IM_Y*state->params.dims.x + IM_BASE_X;
        int copyLen = blockDim.x;

        if(blockIdx.x == 0){
            // We're on the left side, so we have to handle edge conditions
            dstPtr += 1;
            srcPtr += 1;
            copyLen -= 1;
        }else if(blockIdx.x == gridDim.x - 1){
            // We're on the right side, so we have to handle edge conditions
            copyLen -= 1;
        }

        // Perform copy operation
        // memcpy(dstPtr, srcPtr, copyLen*sizeof(FluidVoxel_t));
        _parl_copy(dstPtr, srcPtr, copyLen*sizeof(FluidVoxel_t), blockDim.x, threadIdx.x);
    }

    // Handle edge conditions while we wait for the row to copy
    FluidVoxel_t edgeCondition;

    _setVoxel(
        &edgeCondition, 
        state->params.boundary_velocity.x,
        state->params.boundary_velocity.y,
        1.0f
    );

    if(IM_X >= state->params.dims.x || IM_X < 0 || IM_Y >= state->params.dims.y || IM_Y < 0){
        rowPtr[threadIdx.x] = edgeCondition;
    }

    __syncwarp();
}

static __device__ void _writeBackRow(SimState_t* state, FluidVoxel_t* rowPtr){
    if(IM_Y >= 0 && IM_Y < state->params.dims.y){
        FluidVoxel_t* srcPtr = rowPtr;
        FluidVoxel_t* dstPtr = state->voxels + IM_Y*state->params.dims.x + IM_BASE_X;
        int copyLen = blockDim.x;

        if(blockIdx.x == 0){
            // We're on the left side, so we have to handle edge conditions
            dstPtr += 1;
            srcPtr += 1;
            copyLen -= 1;
        }else if(blockIdx.x == gridDim.x - 1){
            // We're on the right side, so we have to handle edge conditions
            copyLen -= 1;
        }

        // Perform copy operation
        if(threadIdx.x > 0 && threadIdx.x < blockDim.x - 1 && threadIdx.y > 0 && threadIdx.y < blockDim.y - 1){
            // memcpy(dstPtr, srcPtr, copyLen*sizeof(FluidVoxel_t));
            _parl_copy(dstPtr, srcPtr, copyLen*sizeof(FluidVoxel_t), blockDim.x - 2, threadIdx.x - 1);
        }
        __syncwarp();
    }
}

static __device__ void _tileCpy_g2s(SimState_t* state, FluidVoxel_t* tile){
    FluidVoxel_t edgeCondition;

    _setVoxel(
        &edgeCondition, 
        state->params.boundary_velocity.x,
        state->params.boundary_velocity.y,
        1.0f
    );

    if(IM_X >= 0 && IM_X < state->params.dims.x && IM_Y >= 0 && IM_Y < state->params.dims.y){
        tile[MY_TILE_OFFS] = state->voxels[IM_OFFS(state)];
    }else{
        tile[MY_TILE_OFFS] = edgeCondition;
    }
    __syncwarp();
}

static __device__ void _tileCpy_s2g(SimState_t* state, FluidVoxel_t* tile){
    FluidVoxel_t edgeCondition;

    _setVoxel(
        &edgeCondition, 
        state->params.boundary_velocity.x,
        state->params.boundary_velocity.y,
        1.0f
    );

    if(IM_X >= 0 && IM_X < state->params.dims.x && IM_Y >= 0 && IM_Y < state->params.dims.y){
        if(!(threadIdx.x == 0 || threadIdx.x == blockDim.x - 1 || threadIdx.y == 0 || threadIdx.y == blockDim.y - 1)){
            state->voxels[IM_OFFS(state)] = tile[MY_TILE_OFFS];
        }
    }
    __syncwarp();
}

static __device__ void _run(SimState_t* state){
    __shared__ FluidVoxel_t oldVA[TILE_WIDTH*TILE_HEIGHT];
    __shared__ FluidVoxel_t newVA[TILE_WIDTH*TILE_HEIGHT];

    // if(threadIdx.y == 0){
    //     // oldVA[0].lattice_vectors.named.north = 1.0f;
    //     memcpy(&oldVA, state->voxels, 5*sizeof(FluidVoxel_t));
    //     // _parl_copy(oldVA, state->voxels, 5*sizeof(FluidVoxel_t)/sizeof(int32_t), 32, threadIdx.x);
    // }

    __syncthreads();

    FluidVoxel_t* tileRowPtr = &(oldVA[threadIdx.y*TILE_WIDTH]);
    FluidVoxel_t* tileRowPtr_new = &(newVA[threadIdx.y*TILE_WIDTH]);

    // _populateRow(state, tileRowPtr);              // 1: Load data from global memory into oldVA
    _tileCpy_g2s(state, &(oldVA[0]));
    __syncthreads();
    // if(state->frame == 32){
    //     oldVA[MY_TILE_OFFS].velocity.x = (float)threadIdx.x;
    //     oldVA[MY_TILE_OFFS].velocity.y = (float)threadIdx.y;
    // }
    // __syncthreads();
    // No __syncthreads necessary here - see https://github.com/johnMamish/nuce468-w22-fluidsim/issues/3#issuecomment-1062043483
    _collide(state, IM_X, IM_Y, oldVA);           // 2: Collide step
    __syncthreads();                              // 3: Synchronize
    _stream(state, IM_X, IM_Y, oldVA, newVA);     // 4: Stream step
    __syncthreads();                              // 5: Synchronize
    // _barrierBounceBack(state, IM_X, IM_Y, newVA); // 6: Barrier Bounce-Back step
    // __syncthreads();                              // 7: Synchronize
    // // _writeBackRow(state, tileRowPtr_new);                  // 8: Write newVA back to global memory
    _tileCpy_s2g(state, &(newVA[0]));

    if(IM_X == 0 && IM_Y == 0){
        ++(state->frame);
    }
}

__global__ void Opt2Kernel_Full(KERNEL_PARAMS){
    // Run all steps at once
    _run(state);
}

