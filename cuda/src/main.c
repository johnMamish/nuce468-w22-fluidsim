#include <argp.h>
#include <stdio.h>

#include <io/error.h>
#include <io/handle_args.h>
#include <sim/barrier.h>
#include <sim/kernel.h>
#include <sim/sim.h>

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


    // Draw barriers
    printf("%d barriers\n", cArgs.barrier_count);
    for(int i = 0; i < cArgs.barrier_count; ++i){
        SimBarrier_t b = cArgs.barriers[i];
        switch(b.type){
            case SBT_CIRCLE:{
                fseChk_noexit(createBarrier_circle(&state, b.circle.c, b.circle.r), "Error drawing CLI shape #%d (circle)\n", i+1);
            break;}
            case SBT_LINE:{
                fseChk_noexit(createBarrier_line(&state, b.line.p1, b.line.p2), "Error drawing CLI shape #%d (line)\n", i+1);
            break;}
            default:
                fseFail_noexit(FSE_UNKNOWN, "Tried to draw unknown barrier type %d for CLI shape #%d\n", b.type, i+1);
            break;
        }
    }

    // FloatPoint_t triPoints[] = {{50.0,100.0}, {80.0, 130.0}, {80.0, 70.0}};
    // fseChk_noexit(createBarrier_line(&state, triPoints[0], triPoints[1]), "Error drawing triangle");
    // fseChk_noexit(createBarrier_line(&state, triPoints[1], triPoints[2]), "Error drawing triangle");
    // fseChk_noexit(createBarrier_line(&state, triPoints[0], triPoints[2]), "Error drawing triangle");

    // Note that these changes will not take effect until we call this function
    fseChk(syncSimStateToDevice(&state), "Failed to send state to device");

    float avgTime = 0.0f;

    KernelSet_t k = FluidsimKernels.opt1;
    for(int i = 0; i < cArgs.frames; ++i){
        float f;
        fseChk(doFrame(k, &state, &f), "Failure calculating frame %d", i);
        fseChk(syncSimStateToHost(&state), "Failure synchronizing state to host on frame %d", i);
        if(i % cArgs.output_every == 0){
            fseChk(
                writeLogFrame(simFile, &state),
                "Failure writing frame %d to log file", i
            );
        }
        avgTime += f;
    }

    // Note that we only attempt to close the log file if it is not STDERR.
    if(cArgs.output_file){
        fseAssertEq(fclose(simFile), 0, "Could not close output file.");
    }

    avgTime = avgTime / cArgs.frames;

    printf("Kernel execution time (average): %.3f us\n", avgTime*1000);

    // print_args(&cArgs);
}