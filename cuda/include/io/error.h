#pragma once

#include <stdlib.h>

typedef enum FluidsimError{
    FSE_OK = 0,
    FSE_BAD_SIM_STATE = 1,
    FSE_CUDA_ERROR = 2,

    FSE_ASSERTION_FAILURE = -1
} FluidsimError_t;

#define fseAssert(b) {if(!(b)) fseFail(FSE_ASSERTION_FAILURE);}
#define fseAssertEq(a, b) {fseAssert((a) == (b));}
#define fseChk(err) if(err) {fseFail(err);}
#define fseFail(err) {fluidsim_error_exit(err, "in %s (%s line %d)!\n", __func__, __FILE__, __LINE__);}

const char* fluidsim_err_to_str(const FluidsimError_t err);

void fluidsim_error_exit(FluidsimError_t err, const char* format, ...);