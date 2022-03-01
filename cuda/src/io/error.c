#include <io/error.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>

const char* fluidsim_err_to_str(const FluidsimError_t err){
    switch(err){
        case FSE_OK:
            return "FSE_OK";
        case FSE_BAD_SIM_STATE:
            return "FSE_BAD_SIM_STATE";
        case FSE_CUDA_ERROR:
            return "FSE_CUDA_ERROR";
        case FSE_ASSERTION_FAILURE:
            return "FSE_ASSERTION_FAILURE";
    }

    return "UNKNOWN";
}

void fluidsim_error_exit(FluidsimError_t err, const char* format, ...){
    va_list printf_keys;
    va_start(printf_keys, format);

    fprintf(stderr, "ERROR %d (%s): ", (int)err, fluidsim_err_to_str(err));
    vfprintf(stderr, format, printf_keys);

    va_end(printf_keys);

    exit((int)err);
}