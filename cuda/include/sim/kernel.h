#pragma once

#include <sim/sim_struct.h>

#ifndef __global__
    #define __global__
#endif

// __global__ void NaiveKernel(KERNEL_PARAMS);
__global__ void InitializerKernel(SimState_t* state);

__global__ void NaiveKernel_C(KERNEL_PARAMS);
__global__ void NaiveKernel_X(KERNEL_PARAMS);
__global__ void NaiveKernel_S(KERNEL_PARAMS);
__global__ void NaiveKernel_B(KERNEL_PARAMS);

const KernelList_t FluidsimKernels = {
    .naive = NULL
};
