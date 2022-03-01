#include <argp.h>
#include <stdio.h>

#include <io/error.h>
#include <io/handle_args.h>
#include <sim/sim.h>
#include <sim/kernel.h>

// Configuration parameters for argp
const char* argp_program_version = "fluidsim 0.1-GPU";
const char* argp_program_bug_address 
    = "https://github.com/johnMamish/nuce468-w22-fluidsim/issues/new";

int main(int argc, char** argv){
    CLIArgs_t cArgs;

    int res = get_args(&cArgs, argc, argv);
    if(res != 0) return res;

    SimState_t state;
    fseChk(
        initSimState(cArgs.sim, &state),
        "Failed to initialize simulation state."
    );

    FILE* simFile = cArgs.output_file 
        ? fopen(cArgs.output_file, cArgs.append ? "ab" : "wb") 
        : stdout;
    fseAssert(simFile, "Could not open output file for writing.");
    fseChk(initLogFile(simFile, cArgs.sim), "Failed to initialize log file.");
    

    Kernel_t k = FluidsimKernels.naive;
    for(int i = 0; i < cArgs.frames; ++i){
        fseChk(doFrame(k, &state), "Failure calculating frame %d", i);
        fseChk(
            writeLogFrame(simFile, &state),
            "Failure writing frame %d to log file", i
        );
    }

    // Note that we only attempt to close the log file if it is not STDERR.
    if(cArgs.output_file){
        fseAssertEq(fclose(simFile), 0, "Could not close output file.");
    }

    print_args(&cArgs);
}