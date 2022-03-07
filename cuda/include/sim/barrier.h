#pragma once
#include <io/error.h>
#include <sim/sim_struct.h>

/**
 * @brief Mark a pixel in the simulation as a barrier
 * 
 * @param state The simulation state to modify
 * @param p The coordinates of the pixel to set as a barrier
 */
static inline FluidsimError_t createBarrier_point(SimState_t* state, IntPoint_t p){
    if(
        p.x >= state->params.dims.x || p.y >= state->params.dims.y
        || p.x < 0 || p.y < 0
    ){
        return FSE_OUT_OF_BOUNDS;
    }
    int voxelAddr = p.y*state->params.dims.x + p.x;
    state->voxels[voxelAddr].is_barrier = true;
    return FSE_OK;
}

/**
 * @brief Marks all pixels that lie along a line as barrier pixels
 * 
 * @param state Simulation state to modify
 * @param p0 One of the points on the line
 * @param p1 The other point on the line
 */
FluidsimError_t createBarrier_line(SimState_t* state, FloatPoint_t p0, FloatPoint_t p1);