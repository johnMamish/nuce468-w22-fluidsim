#ifndef _SIMULATION_H
#define _SIMULATION_H

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>

/**
 * Holds parameters that are statically determined for each simulation. The application code
 * can just allocate a single simulation parameters struct and fill it out manually.
 *
 * Ideally, these would be read from a config file.
 */
typedef struct simulation_parameters {
    // Width and height of the simulation in fluid voxels
    int xdim;
    int ydim;

    // boundary condition is the same around the whole perimeter.
    // TODO: maybe it would be nice to make this more flexible.
    float boundary_velocity[2];

    float viscosity;
} simulation_parameters_t;

#define LATTICE_DIR_0 0
#define LATTICE_DIR_N 1
#define LATTICE_DIR_S 2
#define LATTICE_DIR_E 3
#define LATTICE_DIR_W 4
#define LATTICE_DIR_NE 5
#define LATTICE_DIR_SE 6
#define LATTICE_DIR_NW 7
#define LATTICE_DIR_SW 8
#define LATTICE_DIR_NUM_DIRS 9

/**
 *
 */
typedef struct fluid_voxel {
    // microscopic densities along each lattice direction
    // Storage order is determined by LATTICE_DIR #defines
    float lattice_vectors[LATTICE_DIR_NUM_DIRS];

    // macroscopic density and velocity
    float rho;
    float u[2];

    float curl;
} fluid_voxel_t;

/**
 * Barrier exists entirely
 */
typedef struct barrier {
    int xdim, ydim;

    // position and velocity
    float x, y;
    float u[2];
    float mass;

    // anchor and spring constant (in 2D)
    float anchor[2];
    float k[2];

    bool* occupancy;
} barrier_t;

/**
 *
 */
typedef struct simulation_state {
    const simulation_parameters_t* params;

    // number of voxels is implied by this->params->xdim and this->params->ydim
    fluid_voxel_t* voxels;

    barrier_t* barriers;
    int num_barriers;
} simulation_state_t;

/**
 * Returns a newly constructed simulation state.
 *
 * @param[in]     params      Parameters that should be used for this simulation. All fields
 *                            must contain valid values before simulation_state_create() is called.
 */
simulation_state_t* simulation_state_create(const simulation_parameters_t* params);

/**
 * Steps the given simulation state by one simulation step.
 */
void step_simulation_state(simulation_state_t* ss);

bool simulation_state_is_stable(simulation_state_t* ss);

int simulation_state_initialize_log_file(FILE* f, simulation_state_t* ss);

/**
 * Writes the current state of the simulation to the log file
 *
 * @param[in,out] f       file pointer to write to; should be opened with "wb"
 * @param[in]     ss      simulation state to write
 *
 * Returns 0 on success, nonzero on failure.
 */
int simulation_state_append_sim_frame_to_log_file(FILE* f, simulation_state_t* ss);
void set_boundary_conditions(simulation_state_t* ss);

barrier_t* barrier_create_manual(int xdim, int ydim, int x, int y, float kx, float ky, float mass, const bool* occupancy);
barrier_t* barrier_create_lines(const int* points, int npoints, float kx, float ky, float mass);
barrier_t* barrier_create_rectangle(int width, int height);
barrier_t* barrier_create_circle(int width, int height);

void barrier_destroy(barrier_t* barr);

void simulation_state_add_barrier(simulation_state_t* ss, const barrier_t* barr);

#endif
