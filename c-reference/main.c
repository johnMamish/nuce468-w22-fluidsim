/**
 * TODOs:
 *   - config file for simulation parameters
 *   - input file for barriers
 *   - barrier support in simulation
 */

/**
 * outputs a binary file consisting of concatenated
 *
 * the first 16 bytes of the file are the width and height of each frame
 *
 * TODO: might be nice to change this to a JSON type file instead.
 */

#include <stdio.h>

#include "simulation.h"

const static int NUM_STEPS = 2000;

int main(int argc, char** argv)
{
    if (argc != 2) {
        fprintf(stderr, "expected 1 argument\r\n");
        return -1;
    }

    FILE* sim_file = fopen(argv[1], "wb");
    if (sim_file == NULL) {
        fprintf(stderr, "failed to open file %s\r\n", argv[1]);
        return -1;
    }

    const int H = 100; const int W = 200;
    static simulation_parameters_t sim_params;
    sim_params.xdim = W;
    sim_params.ydim = H;
    sim_params.boundary_velocity[0] = 0.1; sim_params.boundary_velocity[1] = 0;
    sim_params.viscosity = 0.005;

    simulation_state_t* sim = simulation_state_create(&sim_params);

    simulation_state_initialize_log_file(sim_file, sim);
    simulation_state_append_sim_frame_to_log_file(sim_file, sim);
    for (int i = 0; i < NUM_STEPS; i++) {
        printf("step %05i\n", i);
        set_boundary_conditions(sim);
        for (int j = 0; j < 50; j++) {
            step_simulation_state(sim);
        }

        simulation_state_append_sim_frame_to_log_file(sim_file, sim);
    }
    printf("\r\n");
    fclose(sim_file);

    return 0;
}
