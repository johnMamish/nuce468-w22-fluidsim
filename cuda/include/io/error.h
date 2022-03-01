#pragma once

#include <stdlib.h>

typedef enum FluidsimError{
    FSE_OK = 0,
    FSE_BAD_SIM_STATE = 1,
    FSE_CUDA_ERROR = 2,
    FSE_MALLOC_FAILED = 3,
    FSE_NULL_PTR = 4,
    FSE_FILE_IO_FAILURE = 5,

    FSE_ASSERTION_FAILURE = -1
} FluidsimError_t;

#define fseAssert(b, msg, ...) {if(!(b)) fseFail(FSE_ASSERTION_FAILURE, msg, ##__VA_ARGS__);}
#define fseAssertEq(a, b, msg, ...) {fseAssert((a) == (b), msg, ##__VA_ARGS__);}
#define fseChk(err, msg, ...) if(err) {fseFail(err, msg, ##__VA_ARGS__);}
#define fseFail(err, msg, ...) {_fluidsim_error_exit_macro(err, __LINE__, __FILE__, __func__, msg, ##__VA_ARGS__);}

const char* fluidsim_err_to_str(const FluidsimError_t err);

void _fluidsim_error_exit_macro(FluidsimError_t err, int line, const char* file, const char* func, const char* format, ...);
void fluidsim_error_exit(FluidsimError_t err, const char* format, ...);