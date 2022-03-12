#pragma once

#include <io/error.h>

#define fseCuChk(err) {fseCuChk_inner((err), __func__, __FILE__, __LINE__);}

static inline void fseCuChk_inner(cudaError_t err, const char* func, const char* fil, int line){
    if(err != cudaSuccess){
        fluidsim_error_exit(
            FSE_CUDA_ERROR, 
            "(%s line %d, in %s): CUDA ERROR %d (%s): %s\n", 
            fil, line, func, 
            (int)err, cudaGetErrorName(err), cudaGetErrorString(err)
        );
    }
}