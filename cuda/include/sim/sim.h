#pragma once

#include <stdbool.h>

#include <io/error.h>
#include <sim/sim_struct.h>

FluidsimError_t initSimState(SimParams_t params, SimState_t* toInit);
FluidsimError_t syncSimStateToHost(SimState_t* h_onHost);
FluidsimError_t doFrame(Kernel_t kernel, SimState_t* sim);
FluidsimError_t getCurl(SimState_t* sim, float* curl);
FluidsimError_t initLogFile(FILE* f, SimParams_t p);
FluidsimError_t writeLogFrame(FILE* f, SimState_t* sim);