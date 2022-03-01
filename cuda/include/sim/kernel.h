#pragma once

#include <sim/sim_struct.h>

#ifndef __global__
    #define __global__
#endif

__global__ void NaiveKernel(KERNEL_PARAMS);

const KernelList_t FluidsimKernels = {
    .naive = NaiveKernel
};
