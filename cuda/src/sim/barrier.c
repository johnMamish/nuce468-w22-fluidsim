#include <io/error.h>
#include <sim/barrier.h>
#include <sim/sim_struct.h>

#include <math.h>

#define PRECISION 8

static inline int reduceAndRound(long long x){
    return (int)((x + (1 << (PRECISION-1))) >> PRECISION);
}

FluidsimError_t createBarrier_line(SimState_t* state, FloatPoint_t p0, FloatPoint_t p1){
    float dx = p1.x - p0.x;
    float dy = p1.y - p0.y;

    float steps = abs(dx) > abs(dy) ? abs(dx) : abs(dy);

    float xInc = dx / steps;
    float yInc = dy / steps;

    float x = p0.x;
    float y = p0.y;

    FluidsimError_t e = FSE_OK;

    for(int i = 0; i <= steps; ++i){
        e = createBarrier_point(
            state, 
            (IntPoint_t){(int)round(x), (int)round(y)}
        );
        x += xInc;
        y += yInc;
    }

    return e;
}