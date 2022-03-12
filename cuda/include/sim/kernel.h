#pragma once

#include <sim/sim_struct.h>

#include <stdbool.h>

#ifndef __global__
    #define __global__
#endif

__global__ void InitializerKernel(SimState_t* state);

// Naive KernelSet
__global__ void NaiveKernel_C(KERNEL_PARAMS);
__global__ void NaiveKernel_X(KERNEL_PARAMS);
__global__ void NaiveKernel_S(KERNEL_PARAMS);
__global__ void NaiveKernel_B(KERNEL_PARAMS);

// Opt1 KernelSet
__global__ void Opt1Kernel_Full(KERNEL_PARAMS);

// Opt2 KernelSet
__global__ void Opt2Kernel_Full(KERNEL_PARAMS);

const KernelList_t FluidsimKernels = {
    .naive = (KernelSet_t){
        .tileWidth = TILE_WIDTH,
        .tileHeight = TILE_WIDTH,
        .tileOverlap = 0,
        .full = NULL,
        .collide = NaiveKernel_C,
        .exchange = NaiveKernel_X,
        .stream = NaiveKernel_S,
        .bounceBack = NaiveKernel_B
    },
    .opt1 = (KernelSet_t){
        .tileWidth = 32,
        .tileHeight = 8,
        .tileOverlap = 1,
        .full = Opt1Kernel_Full,
        .collide = NULL,
        .exchange = NULL,
        .stream = NULL,
        .bounceBack = NULL
    },
    .opt2 = (KernelSet_t){
        .tileWidth = 32,
        .tileHeight = 8,
        .tileOverlap = 1,
        .full = Opt2Kernel_Full,
        .collide = NULL,
        .exchange = NULL,
        .stream = NULL,
        .bounceBack = NULL
    }
};
