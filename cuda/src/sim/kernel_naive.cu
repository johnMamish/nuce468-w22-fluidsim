#include <io/error.h>
#include <sim/sim_struct.h>
#include <sim/kernel.h>

#include <sim/_kernel_util.h>

#include <stdio.h>
#include <stdbool.h>

__device__ void _collide(SimState_t* state, int myX, int myY){
    if(myX < state->params.dims.x && myY < state->params.dims.y){
        const float omega = 1.0 / ((3*state->params.viscosity) + 0.5);

        FluidVoxel_t* myVoxel = _index_voxel(state, myX, myY);
        FluidVoxel_t* myNewVoxel = _index_other_voxel(state, myX, myY);

        const static float ux_dirs[] = {0,-1,-1,-1, 0, 1, 1, 1, 0};
        const static float uy_dirs[] = {1, 1, 0,-1,-1,-1, 0, 1, 0};
        // const static float ux_dirs[] = {0.0, 0.0,  0.0, 1.0, -1.0, 1.0,  1.0, -1.0, -1.0};
        // const static float uy_dirs[] = {0.0, 1.0, -1.0, 0.0,  0.0, 1.0, -1.0,  1.0, -1.0};
        float rho = 0.0, ux = 0.0, uy = 0.0;
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

__device__ void _stream(SimState_t* state, int myX, int myY){
    if(myX < state->params.dims.x && myY < state->params.dims.y){
        FluidVoxel_t* myVoxel = _index_voxel(state, myX, myY);
        FluidVoxel_t* myOldVoxel = _index_other_voxel(state, myX, myY);

        myVoxel->density    = myOldVoxel->density;
        myVoxel->velocity.x = myOldVoxel->velocity.x;
        myVoxel->velocity.y = myOldVoxel->velocity.y;

        // Voxels on the edge of the simulation area will end up copying some of
        // their lattice vectors from "out of bounds" - to handle this, we just
        // calculate a single "edge condition voxel" to represent all voxels
        // that are out-of-bounds, and just copy out-of-bounds vectors from this
        // voxel.
        FluidVoxel_t edgeCondition;
        _setVoxel(
            &edgeCondition, 
            state->params.boundary_velocity.x,
            state->params.boundary_velocity.y,
            1.0f
        );

        // ========= Copy lattice vectors =========
        // --> North
        if(myY == 0){
            myVoxel->lattice_vectors.named.north 
                = edgeCondition.lattice_vectors.named.north;
        }else{
            // Copy from voxel to the south
            myVoxel->lattice_vectors.named.north 
                = _index_other_voxel(state, myX, myY - 1)
                    -> lattice_vectors.named.north;
        }

        // --> South
        if(myY == state->params.dims.y - 1){
            myVoxel->lattice_vectors.named.south 
                = edgeCondition.lattice_vectors.named.south;
        }else{
            // Copy from the voxel to the north
            myVoxel->lattice_vectors.named.south
                = _index_other_voxel(state, myX, myY + 1)
                    -> lattice_vectors.named.south;
        }

        // --> East
        if(myX == state->params.dims.x - 1){
            myVoxel->lattice_vectors.named.east 
                = edgeCondition.lattice_vectors.named.east;
        }else{
            // Copy from the voxel to the west
            myVoxel->lattice_vectors.named.east
                = _index_other_voxel(state, myX + 1, myY)
                    -> lattice_vectors.named.east;
        }

        // --> West
        if(myX == 0){
            myVoxel->lattice_vectors.named.west
                = edgeCondition.lattice_vectors.named.west;
        }else{
            // Copy from the voxel to the east
            myVoxel->lattice_vectors.named.west
                = _index_other_voxel(state, myX - 1, myY)
                    -> lattice_vectors.named.west;
        }

        // --> Northeast
        if(myX == 0 || myY == 0){
            myVoxel->lattice_vectors.named.northeast
                = edgeCondition.lattice_vectors.named.northeast;
        }else{
            // Copy from the voxel to the southwest
            myVoxel->lattice_vectors.named.northeast
                = _index_other_voxel(state, myX + 1, myY - 1)
                    -> lattice_vectors.named.northeast;
        }

        // --> Northwest
        if(myX == state->params.dims.x - 1 || myY == 0){
            myVoxel->lattice_vectors.named.northwest
                = edgeCondition.lattice_vectors.named.northwest;
        }else{
            // Copy from the voxel to the southeast
            myVoxel->lattice_vectors.named.northwest
                = _index_other_voxel(state, myX + 1, myY - 1)
                    -> lattice_vectors.named.northwest;
        }

        // --> Southeast
        if(myX == 0 || myY == state->params.dims.y - 1){
            myVoxel->lattice_vectors.named.southeast
                = edgeCondition.lattice_vectors.named.southeast;
        }else{
            // Copy from the voxel to the northwest
            myVoxel->lattice_vectors.named.southeast
                = _index_other_voxel(state, myX - 1, myY + 1)
                    -> lattice_vectors.named.southeast;
        }

        // --> Southwest
        if(myX == state->params.dims.x - 1 || myY == state->params.dims.y - 1){
            myVoxel->lattice_vectors.named.southwest
                = edgeCondition.lattice_vectors.named.southwest;
        }else{
            // Copy from the voxel to the northeast
            myVoxel->lattice_vectors.named.southwest
                = _index_other_voxel(state, myX + 1, myY + 1)
                    -> lattice_vectors.named.southwest;
        }
    }
}

__device__ void _barrierBounceBack(SimState_t* state, int myX, int myY){
    if(myX < state->params.dims.x && myY < state->params.dims.y){
        FluidVoxel_t* myVoxel = _index_voxel(state, myX, myY);
        if(myVoxel->is_barrier){
            // For every lattice vector except the zero vector
            for(int i = 0; i < LV_IM; ++i){
                // Find the voxel in the opposite direction (i.e. if i == LV_S,
                // get the voxel to the north)
                IntPoint_t voxelDelta = _voxel_delta_in_direction(LV_OPPOSITE_DIR_OF(i));
                int targetX = myX + voxelDelta.x;
                int targetY = myY + voxelDelta.y;

                // Make sure we're not trying to bounce-back out of bounds
                if(
                        targetX < 0 || targetX >= state->params.dims.x
                    || targetY < 0 || targetY >= state->params.dims.y
                ){
                    // Bounce back
                    _index_voxel(state, targetX, targetY)
                        -> lattice_vectors.sequence[LV_OPPOSITE_DIR_OF(i)]
                            = myVoxel->lattice_vectors.sequence[i];
                }
            }
        }
    }
}

__global__ void NaiveKernel(KERNEL_PARAMS){
    // Get pixel ID
    int myX = blockIdx.x * blockDim.x + threadIdx.x;
    int myY = blockIdx.y * blockDim.y + threadIdx.y;

    // Collide
    _collide(state, myX, myY);
    __syncthreads();

    // Swap new/old arrays
    if((myX | myY) == 0) _swap_new_old_voxel_arrays(state);
    __syncthreads();

    // Stream
    _stream(state, myX, myY);
    __syncthreads();

    // Handle barrier bounce-back
    _barrierBounceBack(state, myX, myY);
    __syncthreads();
}

