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
    .config_file = NULL,
    .output_file = NULL,
    .append = false,
    .frames = 750,
    .modified = {
        .config_file = false,
        .output_file = false,
        .append = false,
        .frames = false,
        .dims = false,
        .boundary = false,
        .viscosity = false
    }
};


static char argp_doc[] = "Lattice-Boltzmann 2D Fluid Simulator";
static char argp_args_doc[] = "";

static struct argp_option options[] = {
    {"run-len", 'l', "frames", 0, "Number of frames to run the simulator for."},
    {"config-file", 'c', "cfile", 0, "Location of a json-formatted config file to read simulation parameters from. If not provided, default parameters will be used (except when overridden by CLI options)."},
    {"output", 'o', "ofile", 0, "File to write output data to. By default, this will overwrite any existing file (but see `--append`). If not provided, results will be output to STDOUT."},
    {"append", 'a', 0, 0, "If provided, results will be APPENDED to `output` if the file exists, instead of overwriting the file."},
    {0, 0, 0, OPTION_DOC, "SIMULATION PARAMETERS:"},
    {"dims", 'd', "(xDim, yDim)", 0, "Set the dimensions of the simulation."},
    {"boundary", 'b', "(xVel, yVel)", 0, "Set the 'boundary velocity' vector of the simulation."},
    {"viscosity", 'v', "viscosity", 0, "Set the viscosity of the simulation fluid"},
    {0}
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
