#pragma once

#include <argp.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>

#include "sim/sim.h"

/**
 * @brief Struct used to store CLI arguments for the simulation
 * 
 * @param sim A `SimParams_t` struct to hold simulation parameters. Check
 *            <sim/sim.h> for more information on the `SimParams_t` struct.
 * 
 * @param frames Number of frames to run the simulation for
 * @param config_file A config file to load simulation parameters from
 */
typedef struct CLIArgs {
    SimParams_t sim;
    int frames;
    const char* config_file;
    const char* output_file;
    bool append;
    struct {bool config_file; bool output_file; bool frames; bool append; bool dims; bool boundary; bool viscosity;} modified;
} CLIArgs_t;

/**
 * @brief Default command-line arguments.
 * 
 * This gets deep-copied into a CLIArgs struct when argp is initializing.
 */
const static CLIArgs_t DEFAULT_ARGS = (CLIArgs_t){
    .sim = {
        .dims = {
            .x = 200,
            .y = 200
        },
        .boundary_velocity = {
            .x = 0.1,
            .y = 0
        },
        .viscosity = 0.1
    },
    .frames = 750,
    .config_file = NULL,
    .output_file = NULL,
    .append = false,
    .modified = {
        .config_file = false,
        .output_file = false,
        .frames = false,
        .append = false,
        .dims = false,
        .boundary = false,
        .viscosity = false
    }
};




/**
 * @brief Get command-line arguments
 *
 * @param args Pointer to a struct to hold arguments once they are parsed
 */
error_t get_args(CLIArgs_t* args, int argc, char** argv);

/**
 * @brief Print a human-readable representation of the CLIArgs object `args`
 *
 * @param args CLIArgs object to print
 */
void print_args(CLIArgs_t* args);
