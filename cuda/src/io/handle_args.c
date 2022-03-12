#include <argp.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <io/handle_args.h>

static char argp_doc[] = "Lattice-Boltzmann 2D Fluid Simulator";
static char argp_args_doc[] = "";

static struct argp_option options[] = {
    {"run-len", 'l', "frames", 0, "Number of frames to run the simulator for."},
    {"output-every", 'L', "frames", 0, "Only write a frame to output every n frames."},
    // {"config-file", 'c', "cfile", 0, "Location of a json-formatted config file to read simulation parameters from. If not provided, default parameters will be used (except when overridden by CLI options)."},
    {"output", 'o', "ofile", 0, "File to write output data to. By default, this will overwrite any existing file (but see `--append`). If not provided, results will be output to STDOUT."},
    {"append", 'a', 0, 0, "If provided, results will be APPENDED to `output` if the file exists, instead of overwriting the file."},
    {0, 0, 0, OPTION_DOC, "SIMULATION PARAMETERS:"},
    {"dims", 'd', "(xDim, yDim)", 0, "Set the dimensions of the simulation."},
    {"boundary", 'b', "(xVel, yVel)", 0, "Set the 'boundary velocity' vector of the simulation."},
    {"viscosity", 'v', "viscosity", 0, "Set the viscosity of the simulation fluid"},
    {"barrier", 'B', "(type ...)", 0, "Add a barrier to the simulation. See BARRIER TYPES."},
    {0, 0, 0, OPTION_DOC, "BARRIER TYPES"},
    {"(LINE (x1, y1), (x2, y2))", 0, 0, OPTION_DOC, "A direct line between two points (x1, y1) and (x2, y2)."},
    {"(CIRCLE (xc, yc), r)", 0, 0, OPTION_DOC, "A circle centered at (cx, cy) with radius r"},
    {0}
};

const int NUM_VEC_FORMATS_INT = 4;
const char* VEC_FORMATS_INT[] = {
    "(%i, %i)", // <-- This is the most correct one
    "(%i,%i)",
    "%i, %i",
    "%i,%i"
};

const int NUM_VEC_FORMATS_FLOAT = 4;
const char* VEC_FORMATS_FLOAT[] = {
    "(%f, %f)", // <-- This is the most correct one
    "(%f,%f)",
    "%f, %f",
    "%f,%f"
};

const int NUM_BARRIER_FMTS = 10;
#define BARRIER_STR_LEN 6
const char* BARRIER_FORMATS[] = {
    "(%6s (%i, %i), (%i, %i))",
    "(%6s (%i,%i), (%i,%i))",
    "(%6s %i, %i, %i, %i)",
    "(%6s %i,%i, %i,%i)",
    "(%6s %i,%i,%i,%i)",
    "(%6s (%i, %i), %i)",
    "(%6s (%i,%i), %i)",
    "(%6s %i, %i, %i)",
    "(%6s %i,%i, %i)",
    "(%6s %i,%i,%i)",
};

static inline error_t _parse_int_vector(const char* str, int* x, int* y){
    int res = EOF;

    // Try all vector formats we know of, and see if one of them works.
    for(int i = 0; i < NUM_VEC_FORMATS_INT; ++i){
        res = sscanf(str, VEC_FORMATS_INT[i], x, y);
        if(res == 2) break;
    }

    if(res != 2){
        // Can't figure out how to parse this one. Ignore it.
        return EINVAL;
    }

    // One of the formats worked.
    return 0;
}

static inline error_t _parse_float_vector(const char* str, float* x, float* y){
    int res = EOF;

    // Try all vector formats we know of, and see if one of them works.
    for(int i = 0; i < NUM_VEC_FORMATS_FLOAT; ++i){
        res = sscanf(str, VEC_FORMATS_FLOAT[i], x, y);
        if(res == 2) break;
    }

    if(res != 2){
        // Can't figure out how to parse this one. Ignore it.
        return EINVAL;
    }

    // One of the formats worked.
    return 0;
}

static inline error_t _parse_barrier(const char* str, SimBarrierType_t* bType, int* p1, int* p2, int* p3, int* p4){
    char barrier_type[BARRIER_STR_LEN+1];
    int res = 0;
    for(int i = 0; i < NUM_BARRIER_FMTS; ++i){
        res = sscanf(str, BARRIER_FORMATS[i], &barrier_type[0], p1, p2, p3, p4);
        if(res > 3) break;
    }

    if(res <= 3){
        // Can't figure out how to parse this one. Ignore it.
        return EINVAL;
    }

    if(strcmp(&barrier_type[0], "LINE") == 0){
        *bType = SBT_LINE;
        if(res != 5) return EINVAL; // Check specific number of values returned
    }else if(strcmp(&barrier_type[0], "CIRCLE") == 0){
        *bType = SBT_CIRCLE;
        if(res != 4) return EINVAL; // Check specific number of values returned
    }else{
        *bType = SBT_UNKNOWN;
        return EINVAL; // If the barrier type doesn't match anything, always
                       // reject the cli-arg
    }

    // One of the formats worked.
    return 0;
}

/**
 * @brief Parse a single argument
 * 
 * This function is used by `argp` to handle a single command-line argument
 * 
 * @param key Unique identifier for the argument
 * @param arg Value of the argument
 * @param state State information for the argument parser
 * @return error_t Information on any parsing error which may have occurred
 */
static error_t parse_opt(int key, char* arg, struct argp_state* state){
    // This struct allows us to communicate the parsing information to the
    // main program.
    CLIArgs_t* simArgs = (CLIArgs_t*)(state -> input);

    switch(key){
        case ARGP_KEY_INIT:{
            // Called before parsing any arguments. Use this opportunity to
            // initialize the struct fields
            memcpy(simArgs, &DEFAULT_ARGS, sizeof(CLIArgs_t));
        break;}

        case 'c':{
            simArgs -> config_file = arg;
            simArgs -> modified.config_file = true;
        break;}

        case 'o':{
            simArgs -> output_file = arg;
            simArgs -> modified.output_file = true;
        break;}

        case 'a':{
            simArgs -> append = true;
            simArgs -> modified.append = true;
        break;}

        case 'l':{
            int frames;
            int res = sscanf(arg, "%i", &frames);
            if(res != 1){
                argp_usage(state);
            }

            simArgs -> frames = frames;
            simArgs -> modified.frames = true;
        break;}

        case 'L':{
            int outputEvery;
            int res = sscanf(arg, "%i", &outputEvery);
            if(res != 1){
                argp_usage(state);
            }

            simArgs -> output_every = outputEvery;
            simArgs -> modified.output_every = true;
        break;}

        case 'd':{
            int dimX, dimY;
            int res = _parse_int_vector(arg, &dimX, &dimY);
            if(res != 0){
                argp_usage(state);
                // return res;
            }
            simArgs -> sim.dims.x = dimX;
            simArgs -> sim.dims.y = dimY;
            simArgs -> modified.dims = true;
        break;}

        case 'b':{
            float velX, velY;
            int res = _parse_float_vector(arg, &velX, &velY);
            if(res != 0){
                argp_usage(state);
                // return res;
            }
            simArgs -> sim.boundary_velocity.x = velX;
            simArgs -> sim.boundary_velocity.y = velY;
            simArgs -> modified.boundary = true;
        break;}

        case 'v':{
            float visc = 0.0f;
            int res = sscanf(arg, "%f", & visc);
            if(res != 1){
                argp_usage(state);
                // return EINVAL; // Can't parse the float.
            }
            simArgs -> sim.viscosity = visc;
            simArgs -> modified.viscosity = true;
        break;}

        case 'B':{
            int p1, p2, p3, p4;
            SimBarrierType_t bType;

            int res = _parse_barrier(arg, &bType, &p1, &p2, &p3, &p4);
            if(res != 0 || bType == SBT_UNKNOWN){
                argp_usage(state); // Can't parse the barrier
            }

            SimBarrier_t barrier;
            barrier.type = bType;

            if(bType == SBT_CIRCLE){
                barrier.circle.c = (IntPoint_t){p1, p2};
                barrier.circle.r = p3;
            }else if(bType == SBT_LINE){
                barrier.line.p1 = (IntPoint_t){p1, p2};
                barrier.line.p2 = (IntPoint_t){p3, p4};
            }

            if(simArgs->barriers == NULL){
                simArgs->barriers = (SimBarrier_t*) malloc(MAX_BARRIERS*sizeof(SimBarrier_t));
            }

            simArgs->barriers[simArgs->barrier_count] = barrier;
            ++(simArgs->barrier_count);
        break;}

        case ARGP_KEY_ARG:{
            // This program doesn't take unnamed arguments.
            argp_usage(state);
            // return EINVAL;
        break;}

        default:{
            return ARGP_ERR_UNKNOWN;
        }
    }
    return 0;
}


static struct argp argp_cfg = {options, parse_opt, argp_args_doc, argp_doc};

error_t get_args(CLIArgs_t* args, int argc, char** argv){
    return argp_parse(&argp_cfg, argc, argv, 0, 0, args);
}

#define PA_MOD (" [MODIFIED]")
#define PA_NOMOD (" [DEFAULT]")

#define MODSTR(args, key) (args -> modified.key ? PA_MOD : PA_NOMOD)

void print_args(CLIArgs_t* args){
    fprintf(stderr, "Simulation length: %d frames%s\nOutput Every: %d frames%s\nOutput file: %s%s\n  -> Append: %s%s\nSimulation Parameters:\n  -> Dimensions: (%d, %d)%s\n  -> Boundary Velocity: (%.3f, %.3f)%s\n  -> Fluid Viscosity: %.3f%s\n",
        args -> frames     , MODSTR(args, frames),
        args -> output_every, MODSTR(args, output_every),
        args -> output_file, MODSTR(args, output_file),
        (args -> append ? "yes" : "no"), MODSTR(args, append),
        args -> sim.dims.x, args -> sim.dims.y, MODSTR(args, dims),
        args -> sim.boundary_velocity.x, args -> sim.boundary_velocity.y, MODSTR(args, boundary),
        args -> sim.viscosity, MODSTR(args, viscosity)
    );
}
