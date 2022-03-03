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

#define NUM_LATTICE_VECTORS 9
/**
 * @brief A single voxel (actually a pixel since we're in 2D) of fluid.
 *
 * @param lattice_vectors The lattice vectors, accessible by array or named
 *                        struct
 * @param density The overall density of the fluid at this point. 1D.
 * @param velocity The overall velocity of the fluid at this point. 2D.
 * @param curl The rotational motion, or "curl", of the fluid at this point. 1D.
 * @param is_barrier Whether the voxel represents a barrier that obstructs the
 *                   flow of fluid.
 */
typedef struct FluidVoxel {
    /**
     * @brief Makes the lattice vector array accessible via two schemes. Note
     *        that both subfields of this field refer to the SAME DATA, just by
     *        two different methods of access.
     * 
     * @param named The individual vectors accessible by name
     * @param sequence The vectors as an array, which can easily be iterated
     *                 over
     */
    union{
        struct __attribute__((__packed__)) {
            float zero, north, south, east, west,
                northeast, southeast, northwest, southwest;
        } named;
        float sequence[NUM_LATTICE_VECTORS];
    }lattice_vectors;
    float density;
    FloatPoint_t velocity;
    float curl;
    bool is_barrier;
} FluidVoxel_t;

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