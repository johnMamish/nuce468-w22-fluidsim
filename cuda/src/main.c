#include <argp.h>
#include <stdio.h>

#include "io/handle_args.h"
#include "sim/sim.h"

// Configuration parameters for argp
const char* argp_program_version = "fluidsim 0.1-GPU";
const char* argp_program_bug_address = "https://github.com/johnMamish/nuce468-w22-fluidsim/issues/new";

int main(int argc, char** argv){
    CLIArgs_t cArgs;

    int res = get_args(&cArgs, argc, argv);
    if(res != 0) return res;

    // printf("test 2\n");

    print_args(&cArgs);

    test_fun();
}