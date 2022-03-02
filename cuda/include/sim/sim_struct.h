#pragma once

#include <stdbool.h>

#define TILE_WIDTH 32

/**
 * @brief Simple structs to represent a point in 2D space
 */
typedef struct IntPoint {
    int x, y;
} IntPoint_t;
/**
 * @brief Simple structs to represent a point in 2D space
 */
typedef struct FloatPoint {
    float x, y;
} FloatPoint_t;

/**
 * @brief Statically-defined simulation parameters.
 * 
 * @param dims Dimensions of the simulation (x and y)
 * @param boundary_velocity Fluid velocity conditions at the boundaries
 *                          (x-boundary and y-boundary) of the simulation
 * @param viscosity Viscosity of the fluid
 */
typedef struct SimParams {
    IntPoint_t dims;
    FloatPoint_t boundary_velocity;
    float viscosity;
} SimParams_t;

/**
 * @brief A single voxel (actually a pixel since we're in 2D) of fluid.
 *
 * @param lattice_vectors Density in each individual lattice direction (see
 *                        LatticeDirection_t)
 * @param density The overall density of the fluid at this point. 1D.
 * @param velocity The overall velocity of the fluid at this point. 2D.
 * @param curl The rotational motion, or "curl", of the fluid at this point. 1D.
 * @param is_barrier Whether the voxel represents a barrier that obstructs the
 *                   flow of fluid.
 */
typedef struct FluidVoxel {
    struct {
        float zero, north, south, east, west,
              northeast, southeast, northwest, southwest;
        } lattice_vectors;
    float density;
    FloatPoint_t velocity;
    float curl;
    bool is_barrier;
} FluidVoxel_t;
#define NUM_LATTICE_VECTORS 9

/**
 * @brief Simulation state information
 * 
 * @param frame The current frame of the simulation
 * @param params The parameters of the simulation
 * @param voxels An array of voxels
 * 
 * N.B. - Barrier information is now encoded within the voxels themselves,
 * rather than in a separate array.
 */
typedef struct SimState {
    int frame;
    SimParams_t   params;
    FluidVoxel_t* voxels;
    struct SimState* d_deviceStatePtr;
} SimState_t;

#define KERNEL_PARAMS SimState_t* state

typedef void (*Kernel_t)(KERNEL_PARAMS);

typedef struct KernelList {
    Kernel_t naive;
} KernelList_t;