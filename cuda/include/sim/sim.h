#pragma once

/* ======== Data Structures ======== */

/**
 * @brief Statically-defined simulation parameters.
 * 
 * @param dims Dimensions of the simulation (x and y)
 * @param boundary_velocity Fluid velocity conditions at the boundaries
 *                          (x-boundary and y-boundary) of the simulation
 * @param viscosity Viscosity of the fluid
 */
typedef struct SimParams {
    struct {int x; int y;} dims;
    struct {float x; float y;} boundary_velocity;
    float viscosity;
} SimParams_t;

/**
 * @brief D2Q9 lattice directions
 * 
 * This simulator uses a D2Q9 vector configuration for the lattice-boltzmann
 * simulation - this enum defines human-readable names for the different
 * vectors, in order to ease debugging and increase code readability.
 */
enum LatticeDirection_t {
    LD_0  = 0, // The zero vector
    LD_N  = 1, // Vector in the "north" (up) direction
    LD_S  = 2, // Vector in the "south" (down) direction
    LD_E  = 3, // Vector in the "east" (right) direction
    LD_W  = 4, // Vector in the "west" (left) direction
    LD_NE = 5, // Vector in the "northeast" (up-right) direction
    LD_SE = 6, // Vector in the "southeast" (down-right) direction
    LD_NW = 7, // Vector in the "northwest" (up-left) direction
    LD_SW = 8, // Vector in the "southwest" (down-left) direction
};

/* ======== Function Prototypes ======== */

void test_fun();