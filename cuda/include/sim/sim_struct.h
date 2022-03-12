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

struct SimBarrier_line {
    IntPoint_t p1, p2;
};

struct SimBarrier_circle {
    IntPoint_t c;
    int r;
};

typedef enum {
    SBT_CIRCLE,
    SBT_LINE,
    SBT_UNKNOWN,
} SimBarrierType_t;

typedef struct SimBarrier {
    SimBarrierType_t type;
    struct SimBarrier_line line;
    struct SimBarrier_circle circle;
} SimBarrier_t;

#define NUM_LATTICE_VECTORS 9
#define LV_IP 4
#define LV_IM 8

#define LV_N  0
#define LV_NW 1
#define LV_W  2
#define LV_SW 3
#define LV_S  4
#define LV_SE 5
#define LV_E  6
#define LV_NE 7
#define LV_Z  8
// NOTE: this doesn't work for the zero vector.
#define LV_OPPOSITE_DIR_OF(x) ((x + LV_IP) % LV_IM)
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
            // NOTE: The exact ordering of vectors here is VERY important, as
            // it makes it possible to invert a direction (i.e. north -> south,
            // northeast -> southwest) by using the following calculation:
            // int newDirIndex = (oldDirIndex + LV_IP) % LV_IM;
            float north, northwest, west, southwest,
                  south, southeast, east, northeast,
                  zero;
        } named;
        float sequence[NUM_LATTICE_VECTORS];
    }lattice_vectors;
    float density;
    FloatPoint_t velocity;
    float curl;
    bool is_barrier;
    int _pad; // Pads struct from 14*4 to 15*4, to avoid shared memory conflicts
} FluidVoxel_t;

/**
 * @brief Simulation state information
 *
 * @param frame The current frame of the simulation
 * @param params The parameters of the simulation
 * @param voxels An array of voxels
 * @param d_deviceStatePtr Pointer to the address in which the simulation state
 *                         is stored on the device.
 *
 * @param _d_voxels_old This pointer is only valid on device. It is used
 *                      internally by the Lattice-Boltzmann algorithm, and
 *                      should never be touched by the host.
 *
 * N.B. - Barrier information is now encoded within the voxels themselves,
 * rather than in a separate array.
 */
typedef struct SimState {
    int frame;
    SimParams_t   params;
    FluidVoxel_t* voxels;
    struct SimState* d_deviceStatePtr;
    FluidVoxel_t* _d_voxels_old;
} SimState_t;

#define KERNEL_PARAMS SimState_t* state

typedef void (*Kernel_t)(KERNEL_PARAMS);

typedef struct KernelSet {
    unsigned int tileWidth, tileHeight, tileOverlap;
    Kernel_t full, collide, exchange, stream, bounceBack;
} KernelSet_t;

typedef struct KernelList {
    KernelSet_t naive, opt1, opt2;
} KernelList_t;