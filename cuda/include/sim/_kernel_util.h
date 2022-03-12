#pragma once

#include <io/error.h>
#include <sim/sim_struct.h>
#include <sim/kernel.h>

#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>

#define _VOXEL_IDX(sim,xx,yy)(yy*(sim->params.dims.x) + xx)

/**
 * @brief Gets the (x,y) offsets corresponding to the direction provided
 * 
 * e.g. _voxel_delta_in_direction(LV_NE) returns {1, 1} which, if you added
 * it to a voxel x and y coordinate, would return the voxel northeast of the
 * base voxel (i.e. the one with x+1 and y+1 to the base voxel)
 * 
 * @param lv_dir The lattice-vector direction. use LV_XX defines in sim_struct.h
 * @return {x_offset, y_offset}
 */
static inline __device__ IntPoint_t _voxel_delta_in_direction(int lv_dir){
    switch(lv_dir){
        case LV_N:
            return (IntPoint_t){ 0,  1};
        case LV_S:
            return (IntPoint_t){ 0, -1};
        case LV_E:
            return (IntPoint_t){ 1,  0};
        case LV_W:
            return (IntPoint_t){-1,  0};
        case LV_NE:
            return (IntPoint_t){ 1,  1};
        case LV_SE:
            return (IntPoint_t){ 1, -1};
        case LV_NW:
            return (IntPoint_t){-1,  1};
        case LV_SW:
            return (IntPoint_t){-1, -1};
        case LV_Z:
            return (IntPoint_t){ 0,  0};
        default:
            return (IntPoint_t){ 0,  0};
    }
}

static inline __device__ FluidVoxel_t* _index_voxel(SimState_t* sim, int x, int y)
{
    return &(sim->voxels[_VOXEL_IDX(sim,x,y)]);
}

static inline __device__ FluidVoxel_t* _index_other_voxel(SimState_t* sim, int x, int y)
{
    return &(sim->_d_voxels_old[_VOXEL_IDX(sim,x,y)]);
}

static inline __device__ void _swap_new_old_voxel_arrays(SimState_t* sim){
    FluidVoxel_t* nextOld = sim->voxels;
    sim->voxels = sim->_d_voxels_old;
    sim->_d_voxels_old = nextOld;
}

static inline __device__ void _mod_voxel_lattice_vectors(
    FluidVoxel_t* v,
    float xVel,
    float yVel,
    float density,
    float omega
){
    // Calculate intermediate values for lattice vectors
    float ux_times_3 = 3.0f  * xVel;
    float uy_times_3 = 3.0f  * yVel;
    float ux_squared = xVel * xVel;
    float uy_squared = yVel * yVel;
    float ux_times_uy_times_2 = xVel * yVel * 2.0f;
    float u_squared = ux_squared + uy_squared;
    float u_squared_times_150pct = u_squared * 1.5f;

    float zero_mag = 4.0f/9.0f;
    float nesw_mag = 1.0f/9.0f;
    float crnr_mag = 1.0f/36.0f;

    // Calculate lattice vectors
    v->lattice_vectors.named.zero += omega * (
        zero_mag
        * (1.0f - u_squared_times_150pct) 
        * density 
        - v->lattice_vectors.named.zero);

    v->lattice_vectors.named.east  += omega * (
        nesw_mag
        * (1.0f + ux_times_3 + 4.5f*ux_squared - u_squared_times_150pct) 
        * density 
        - v->lattice_vectors.named.east
    );
    v->lattice_vectors.named.west  += omega * (
        nesw_mag
        *(1.0f - ux_times_3 + 4.5f*ux_squared - u_squared_times_150pct) 
        * density 
        - v->lattice_vectors.named.west
    );
    v->lattice_vectors.named.north += omega * (
        nesw_mag
        *(1.0f + uy_times_3 + 4.5f*uy_squared - u_squared_times_150pct) 
        * density 
        - v->lattice_vectors.named.north
    );
    v->lattice_vectors.named.south += omega * (
        nesw_mag
        *(1.0f - uy_times_3 + 4.5f*uy_squared - u_squared_times_150pct) 
        * density 
        - v->lattice_vectors.named.south
    );

    v->lattice_vectors.named.northeast += omega * (
        crnr_mag
        *(1.0f + ux_times_3 + uy_times_3 
            + 4.5f*(u_squared + ux_times_uy_times_2) 
            - u_squared_times_150pct
        )
        * density 
        - v->lattice_vectors.named.northeast
    );
    v->lattice_vectors.named.southeast += omega * (
        crnr_mag
        *(1.0f + ux_times_3 - uy_times_3 
            + 4.5f*(u_squared - ux_times_uy_times_2) 
            - u_squared_times_150pct)
        * density 
        - v->lattice_vectors.named.southeast
    );
    v->lattice_vectors.named.northwest += omega * (
        crnr_mag
        *(1.0f - ux_times_3 + uy_times_3 
            + 4.5f*(u_squared - ux_times_uy_times_2) 
            - u_squared_times_150pct)
        * density 
        - v->lattice_vectors.named.northwest
    );
    v->lattice_vectors.named.southwest += omega * (
        crnr_mag
        *(1.0f - ux_times_3 - uy_times_3 
            + 4.5f*(u_squared + ux_times_uy_times_2) 
            - u_squared_times_150pct)
        * density 
        - v->lattice_vectors.named.southwest
    );

}

static inline __device__ void _setVoxel(FluidVoxel_t* v, float xVel, float yVel, float density){
    // Set lattice vectors
    for(int i = 0; i < NUM_LATTICE_VECTORS; ++i) v->lattice_vectors.sequence[i] = 0.0f;
    _mod_voxel_lattice_vectors(v, xVel, yVel, density, 1.0f);

    // Copy over other values
    v->density = density;
    v->velocity.x = xVel;
    v->velocity.y = yVel;
}