#include <io/error.h>
#include <sim/sim_struct.h>
#include <sim/kernel.h>

#include <sim/_kernel_util.h>

#include <stdio.h>
#include <stdbool.h>

__global__ void InitializerKernel(SimState_t* state){
    // Initialize simulation state in parallel
    int myX = blockIdx.x * blockDim.x + threadIdx.x;
    int myY = blockIdx.y * blockDim.y + threadIdx.y;

    float imSize = state->params.dims.x * state->params.dims.y;

    if(myX < state->params.dims.x && myY < state->params.dims.y){
        FluidVoxel_t* v = _index_voxel(state, myX, myY);
        FluidVoxel_t* v2 = _index_other_voxel(state, myX, myY);
        _setVoxel(
            v, 
            state->params.boundary_velocity.x,
            state->params.boundary_velocity.y,
            1.0
        );
        _setVoxel(
            v2, 
            state->params.boundary_velocity.x,
            state->params.boundary_velocity.y,
            1.0
        );
        // _setVoxel(
        //     v, 
        //     state->params.boundary_velocity.x,
        //     state->params.boundary_velocity.y,
        //     1.0
        // );
        // _setVoxel(
        //     v2, 
        //     state->params.boundary_velocity.x,
        //     state->params.boundary_velocity.y,
        //     1.0
        // );
        v->curl = 0.0;
        v2->curl = 0.0;
        v->is_barrier = false;
        v2->is_barrier = false;
    }
}