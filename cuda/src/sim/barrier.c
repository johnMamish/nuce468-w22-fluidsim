#include <io/error.h>
#include <sim/barrier.h>
#include <sim/sim_struct.h>

#include <math.h>

#define PRECISION 8

static inline int reduceAndRound(long long x){
    return (int)((x + (1 << (PRECISION-1))) >> PRECISION);
}

FluidsimError_t createBarrier_line(SimState_t* state, IntPoint_t p0, IntPoint_t p1){
    FloatPoint_t p0_f = {(float)p0.x, (float)p0.y};
    FloatPoint_t p1_f = {(float)p1.x, (float)p1.y};

    float dx = p1_f.x - p0_f.x;
    float dy = p1_f.y - p0_f.y;

    float steps = abs(dx) > abs(dy) ? abs(dx) : abs(dy);

    float xInc = dx / steps;
    float yInc = dy / steps;

    float x = p0_f.x;
    float y = p0_f.y;

    int e = FSE_OK;

    for(int i = 0; i <= steps; ++i){
        e |= createBarrier_point(
            state, 
            (IntPoint_t){(int)round(x), (int)round(y)}
        );
        x += xInc;
        y += yInc;
    }

    return (FluidsimError_t)e;
}

// Geometry breaks my brain
// Thanks to https://www.geeksforgeeks.org/mid-point-circle-drawing-algorithm/
// for the circle-drawing algorithm
static inline FluidsimError_t _createReflectedPoints(
    SimState_t* state, 
    IntPoint_t c, 
    IntPoint_t p
){
    int e = FSE_OK;
    e |= createBarrier_point(state, (IntPoint_t){c.x + p.x, c.y + p.y});
    e |= createBarrier_point(state, (IntPoint_t){c.x - p.x, c.y + p.y});
    e |= createBarrier_point(state, (IntPoint_t){c.x + p.x, c.y - p.y});
    e |= createBarrier_point(state, (IntPoint_t){c.x - p.x, c.y - p.y});

    e |= createBarrier_point(state, (IntPoint_t){c.x + p.y, c.y + p.x});
    e |= createBarrier_point(state, (IntPoint_t){c.x - p.y, c.y + p.x});
    e |= createBarrier_point(state, (IntPoint_t){c.x + p.y, c.y - p.x});
    e |= createBarrier_point(state, (IntPoint_t){c.x - p.y, c.y - p.x});

    return (FluidsimError_t)e;
}

FluidsimError_t createBarrier_circle(SimState_t* state, IntPoint_t c, int r){

    int x = r;
    int y = 0;

    int e = FSE_OK;

    // Initial points
    e |= createBarrier_point(state, (IntPoint_t){c.x + r, c.y});
    e |= createBarrier_point(state, (IntPoint_t){c.x - r, c.y});
    e |= createBarrier_point(state, (IntPoint_t){c.x, c.y + r});
    e |= createBarrier_point(state, (IntPoint_t){c.x, c.y - r});

    int p = 1 - r;
    while(x > y){
        ++y;

        if(p <= 0){
            p += 2*y + 1;
        }else{
            --x;
            p += 2*y - 2*x + 1;
        }

        e |= _createReflectedPoints(state, c, (IntPoint_t){x,y});
    }

    return (FluidsimError_t)e;
}