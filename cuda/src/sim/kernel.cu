#include <io/error.h>
#include <sim/sim_struct.h>
#include <sim/kernel.h>

#include <stdio.h>
#include <stdbool.h>

static inline __device__ FluidVoxel_t* index_voxel(SimState_t* sim, int x, int y)
{
    return &(sim->voxels[y * (sim->params.dims.x) + x]);
}

__device__ void _setVoxel(FluidVoxel_t* v, float xVel, float yVel, float density){
    // Calculate intermediate values for lattice vectors
    float ux_times_3 = 3.0  * xVel;
    float uy_times_3 = 3.0  * yVel;
    float ux_squared = xVel * xVel;
    float uy_squared = yVel * yVel;
    float ux_times_uy_times_2 = xVel * yVel * 2.0;
    float u_squared = ux_squared + uy_squared;
    float u_squared_times_150pct = u_squared * 1.5;

    float zero_mag = 4.0/9.0;
    float nesw_mag = 1.0/9.0;
    float crnr_mag = 1.0/36.0;

    // Calculate lattice vectors
    v->lattice_vectors.zero = zero_mag*(1 - u_squared_times_150pct) * density;

    v->lattice_vectors.east  = nesw_mag*(1 + ux_times_3 + 4.5*ux_squared - u_squared_times_150pct) * density;
    v->lattice_vectors.west  = nesw_mag*(1 - ux_times_3 + 4.5*ux_squared - u_squared_times_150pct) * density;
    v->lattice_vectors.north = nesw_mag*(1 + uy_times_3 + 4.5*uy_squared - u_squared_times_150pct) * density;
    v->lattice_vectors.south = nesw_mag*(1 - uy_times_3 + 4.5*uy_squared - u_squared_times_150pct) * density;

    v->lattice_vectors.northeast = crnr_mag*(1 + ux_times_3 + uy_times_3 + 4.5*(u_squared + ux_times_uy_times_2) - u_squared_times_150pct) * density;
    v->lattice_vectors.southeast = crnr_mag*(1 + ux_times_3 - uy_times_3 + 4.5*(u_squared - ux_times_uy_times_2) - u_squared_times_150pct) * density;
    v->lattice_vectors.northwest = crnr_mag*(1 - ux_times_3 + uy_times_3 + 4.5*(u_squared - ux_times_uy_times_2) - u_squared_times_150pct) * density;
    v->lattice_vectors.southwest = crnr_mag*(1 - ux_times_3 - uy_times_3 + 4.5*(u_squared + ux_times_uy_times_2) - u_squared_times_150pct) * density;

    // Copy over other values
    v->density = density;
    v->velocity.x = xVel;
    v->velocity.y = yVel;
}

__global__ void InitializerKernel(SimState_t* state){
    // Initialize simulation state in parallel
    int myX = blockIdx.x * blockDim.x + threadIdx.x;
    int myY = blockIdx.y * blockDim.y + threadIdx.y;

    float imSize = state->params.dims.x * state->params.dims.y;

    if(myX < state->params.dims.x && myY < state->params.dims.y){
        FluidVoxel_t* v = index_voxel(state, myX, myY);
        _setVoxel(
            v, 
            (float)myX,
            (float)myY*-1.0f,
            1.0
        );
        // _setVoxel(
        //     v, 
        //     state->params.boundary_velocity.x,
        //     state->params.boundary_velocity.y,
        //     1.0
        // );
        v->curl = 0.0;
    }
}

__global__ void NaiveKernel(KERNEL_PARAMS){
    // Right now, the kernel just displays a test pattern
    int myX = blockIdx.x * blockDim.x + threadIdx.x;
    int myY = blockIdx.y * blockDim.y + threadIdx.y;

    if(myX < state->params.dims.x && myY < state->params.dims.y){
        FluidVoxel_t* v = index_voxel(state, myX, myY);
        if(myY*state->params.dims.x + myX <= 5*state->frame){
            _setVoxel(
                v, 
                (float)myY,
                (float)myX,
                1.0
            );
        }
    }
    int frame = state->frame;
    __syncthreads();
    state->frame = frame + 1;
}

