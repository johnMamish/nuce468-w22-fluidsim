/**
 * Code is adapted from the javascript found at https://physics.weber.edu/schroeder/fluids/
 * originally by Dan Schroeder
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "simulation.h"

static fluid_voxel_t* index_voxel(simulation_state_t* ss, int x, int y)
{
    return &ss->voxels[y * (ss->params->xdim) + x];
}

/**
 * Initializes a voxel's lattice vectors, density, and static.
 */
static void set_voxel(fluid_voxel_t* v, const float* u_0, const float rho_0)
{
    float ux3 = 3 * u_0[0];
    float uy3 = 3 * u_0[1];
    float ux2 = u_0[0] * u_0[0];
    float uy2 = u_0[1] * u_0[1];
    float uxuy2 = 2 * u_0[0] * u_0[1];
    float u2 = ux2 + uy2;
    float u215 = 1.5 * u2;
    v->lattice_vectors[LATTICE_DIR_0] = (4.0 / 9.0) * (1                              - u215);
    v->lattice_vectors[LATTICE_DIR_N] = (1.0 / 9.0) * (1 + ux3       + 4.5*ux2        - u215);
    v->lattice_vectors[LATTICE_DIR_S] = (1.0 / 9.0) * (1 - ux3       + 4.5*ux2        - u215);
    v->lattice_vectors[LATTICE_DIR_E] = (1.0 / 9.0) * (1 + uy3       + 4.5*uy2        - u215);
    v->lattice_vectors[LATTICE_DIR_W] = (1.0 / 9.0) * (1 - uy3       + 4.5*uy2        - u215);
    v->lattice_vectors[LATTICE_DIR_NE] =(1.0 / 36.0) * (1 + ux3 + uy3 + 4.5*(u2+uxuy2) - u215);
    v->lattice_vectors[LATTICE_DIR_SE] =(1.0 / 36.0) * (1 + ux3 - uy3 + 4.5*(u2-uxuy2) - u215);
    v->lattice_vectors[LATTICE_DIR_NW] =(1.0 / 36.0) * (1 - ux3 + uy3 + 4.5*(u2-uxuy2) - u215);
    v->lattice_vectors[LATTICE_DIR_SW] =(1.0 / 36.0) * (1 - ux3 - uy3 + 4.5*(u2+uxuy2) - u215);

    for (int i = 0; i < LATTICE_DIR_NUM_DIRS; i++)
        v->lattice_vectors[i] *= rho_0;

    v->rho = rho_0;
    v->u[0] = u_0[0];
    v->u[1] = u_0[1];
}

void dump_simulation_state(simulation_state_t* ss) {
    for (int y = 0; y < ss->params->ydim; y++) {
        for (int x = 0; x < ss->params->xdim; x++) {
            printf("%2.1f ", (double)ss->voxels[(y * ss->params->xdim) + x].rho);
        }
        printf("\n");
    }
    printf("\n");
}

void set_boundary_conditions(simulation_state_t* ss)
{
    fluid_voxel_t* v;
    for (int x = 0; x < ss->params->xdim; x++) {
        v = index_voxel(ss, x, 0);
        set_voxel(v, ss->params->boundary_velocity, 1.0);

        v = index_voxel(ss, x, (ss->params->ydim - 1));
        set_voxel(v, ss->params->boundary_velocity, 1.0);
    }

    for (int y = 0; y < ss->params->ydim; y++) {
        v = index_voxel(ss, 0, y);
        set_voxel(v, ss->params->boundary_velocity, 1.0);

        v = index_voxel(ss, (ss->params->xdim - 1), y);
        set_voxel(v, ss->params->boundary_velocity, 1.0);
    }
}

simulation_state_t* simulation_state_create(const simulation_parameters_t* params)
{
    simulation_state_t* ss = calloc(1, sizeof(simulation_state_t));

    ss->params = params;
    ss->voxels = calloc(params->xdim * params->ydim, sizeof(fluid_voxel_t));
    ss->barriers = NULL;
    ss->num_barriers = 0;

    // initialize the whole thing to the boundary condition
    for (int y = 0; y < params->ydim; y++) {
        for (int x = 0; x < params->xdim; x++) {
            fluid_voxel_t* v = index_voxel(ss, x, y);
            set_voxel(v, ss->params->boundary_velocity, 1.0);
            v->curl = 0.0;
        }
    }

    return ss;
}

static void collide(simulation_state_t* ss)
{
    const float omega = 1.0 / ((3 * ss->params->viscosity) + 0.5);      // reciprocal of relaxation time
    for (int y = 1; y < (ss->params->ydim - 1); y++) {
        for (int x = 1; x < (ss->params->xdim - 1); x++) {
            fluid_voxel_t* v = index_voxel(ss, x, y);

            const static char ux_dirs[] = {0,  0,  0,  1, -1,  1,  1, -1, -1};
            const static char uy_dirs[] = {0,  1, -1,  0,  0,  1, -1,  1, -1};
            float thisrho = 0., thisux = 0., thisuy = 0.;
            for (int i = 0; i < LATTICE_DIR_NUM_DIRS; thisrho += v->lattice_vectors[i++]);
            for (int i = 0; i < LATTICE_DIR_NUM_DIRS; i++)
                thisux += (float)ux_dirs[i] * v->lattice_vectors[i];
            for (int i = 0; i < LATTICE_DIR_NUM_DIRS; i++)
                thisuy += (float)uy_dirs[i] * v->lattice_vectors[i];

            thisux /= thisrho;
            thisuy /= thisrho;

            v->rho = thisrho;
            v->u[0] = thisux;
            v->u[1] = thisuy;

            float ux3 = 3 * thisux;
            float uy3 = 3 * thisuy;
            float ux2 = thisux * thisux;
            float uy2 = thisuy * thisuy;
            float uxuy2 = 2 * thisux * thisuy;
            float u2 = ux2 + uy2;
            float u215 = 1.5 * u2;
            v->lattice_vectors[LATTICE_DIR_0]  += omega * ((4. / 9.) * thisrho *
                                                           (1                              - u215) -
                                                           v->lattice_vectors[LATTICE_DIR_0]);
            v->lattice_vectors[LATTICE_DIR_E]  += omega * ((1. / 9.) * thisrho *
                                                           (1 + ux3       + 4.5*ux2        - u215) -
                                                           v->lattice_vectors[LATTICE_DIR_E]);
            v->lattice_vectors[LATTICE_DIR_W]  += omega * ((1. / 9.) * thisrho *
                                                           (1 - ux3       + 4.5*ux2        - u215) -
                                                           v->lattice_vectors[LATTICE_DIR_W]);
            v->lattice_vectors[LATTICE_DIR_N]  += omega * ((1. / 9.) * thisrho *
                                                           (1 + uy3       + 4.5*uy2        - u215) -
                                                           v->lattice_vectors[LATTICE_DIR_N]);
            v->lattice_vectors[LATTICE_DIR_S]  += omega * ((1. / 9.) * thisrho *
                                                           (1 - uy3       + 4.5*uy2        - u215) -
                                                           v->lattice_vectors[LATTICE_DIR_S]);
            v->lattice_vectors[LATTICE_DIR_NE] += omega * ((1. / 36.) * thisrho *
                                                           (1 + ux3 + uy3 + 4.5*(u2+uxuy2) - u215) -
                                                           v->lattice_vectors[LATTICE_DIR_NE]);
            v->lattice_vectors[LATTICE_DIR_SE] += omega * ((1. / 36.) * thisrho *
                                                           (1 + ux3 - uy3 + 4.5*(u2-uxuy2) - u215) -
                                                           v->lattice_vectors[LATTICE_DIR_SE]);
            v->lattice_vectors[LATTICE_DIR_NW] += omega * ((1. / 36.) * thisrho *
                                                           (1 - ux3 + uy3 + 4.5*(u2-uxuy2) - u215) -
                                                           v->lattice_vectors[LATTICE_DIR_NW]);
            v->lattice_vectors[LATTICE_DIR_SW] += omega * ((1. / 36.) * thisrho *
                                                           (1 - ux3 - uy3 + 4.5*(u2+uxuy2) - u215) -
                                                           v->lattice_vectors[LATTICE_DIR_SW]);
        }
    }

    // "at right end, copy left-flowing densities from next row to the left"
    // TODO: will this loop be different for different boundary conditions?
    for (int y = 1; y < (ss->params->ydim - 2); y++) {
        fluid_voxel_t* v_left = index_voxel(ss, ss->params->xdim - 2, y);
        fluid_voxel_t* v = index_voxel(ss, ss->params->xdim - 1, y);

        v->lattice_vectors[LATTICE_DIR_W] = v_left->lattice_vectors[LATTICE_DIR_W];
        v->lattice_vectors[LATTICE_DIR_NW] = v_left->lattice_vectors[LATTICE_DIR_NW];
        v->lattice_vectors[LATTICE_DIR_SW] = v_left->lattice_vectors[LATTICE_DIR_SW];
    }
}

static void stream(simulation_state_t* ss)
{
    //float barrierCount = 0, barrierxSum = 0, barrierySum = 0;
    //float barrierFx = 0.0, barrierFy = 0.0;

    // first start in NW corner...
    for (int y = (ss->params->ydim - 2); y > 0; y--) {
        for (int x = 1; (x < ss->params->xdim - 1); x++) {
            fluid_voxel_t* v = index_voxel(ss, x, y);

            v->lattice_vectors[LATTICE_DIR_N] = index_voxel(ss, x, y - 1)->lattice_vectors[LATTICE_DIR_N];
            v->lattice_vectors[LATTICE_DIR_NW] = index_voxel(ss, x + 1, y - 1)->lattice_vectors[LATTICE_DIR_NW];
        }
    }

    // now start in NE corner...
    for (int y = (ss->params->ydim - 2); y > 0; y--) {
        for (int x = (ss->params->xdim - 2); x > 0; x--) {
            fluid_voxel_t* v = index_voxel(ss, x, y);

            v->lattice_vectors[LATTICE_DIR_E] = index_voxel(ss, x - 1, y)->lattice_vectors[LATTICE_DIR_E];
            v->lattice_vectors[LATTICE_DIR_NE] = index_voxel(ss, x - 1, y - 1)->lattice_vectors[LATTICE_DIR_NE];
        }
    }

    // now start in SE corner...
    for (int y = 1; y < (ss->params->ydim - 1); y++) {
        for (int x = (ss->params->xdim - 2); x > 0; x--) {
            fluid_voxel_t* v = index_voxel(ss, x, y);

            v->lattice_vectors[LATTICE_DIR_S] = index_voxel(ss, x, y + 1)->lattice_vectors[LATTICE_DIR_S];
            v->lattice_vectors[LATTICE_DIR_SE] = index_voxel(ss, x - 1, y + 1)->lattice_vectors[LATTICE_DIR_SE];
        }
    }

    // now start in the SW corner...
    for (int y = 1; y < (ss->params->ydim - 1); y++) {
        for (int x=1; x < (ss->params->xdim - 1); x++) {
            fluid_voxel_t* v = index_voxel(ss, x, y);

            v->lattice_vectors[LATTICE_DIR_W] = index_voxel(ss, x + 1, y)->lattice_vectors[LATTICE_DIR_W];
            v->lattice_vectors[LATTICE_DIR_SW] = index_voxel(ss, x + 1, y + 1)->lattice_vectors[LATTICE_DIR_SW];
        }
    }

#if 1
    // Now handle bounce-back from barriers
    for (int y = 1; y < (ss->params->ydim - 1); y++) {
        for (int x = 1; x < (ss->params->xdim - 1); x++) {
            /*if (barrier[x+y*xdim]) { */
            if ((x == 60) && (y >= 25) && (y < 75)) {
                //if (0) {
                fluid_voxel_t* v = index_voxel(ss, x, y);

                index_voxel(ss, x + 1, y)->lattice_vectors[LATTICE_DIR_E] = v->lattice_vectors[LATTICE_DIR_W];
                index_voxel(ss, x - 1, y)->lattice_vectors[LATTICE_DIR_W] = v->lattice_vectors[LATTICE_DIR_E];
                index_voxel(ss, x, y + 1)->lattice_vectors[LATTICE_DIR_N] = v->lattice_vectors[LATTICE_DIR_S];
                index_voxel(ss, x, y - 1)->lattice_vectors[LATTICE_DIR_S] = v->lattice_vectors[LATTICE_DIR_N];

                index_voxel(ss, x + 1, y + 1)->lattice_vectors[LATTICE_DIR_NE] =
                    v->lattice_vectors[LATTICE_DIR_SW];
                index_voxel(ss, x - 1, y + 1)->lattice_vectors[LATTICE_DIR_NW] =
                    v->lattice_vectors[LATTICE_DIR_SE];
                index_voxel(ss, x + 1, y - 1)->lattice_vectors[LATTICE_DIR_SE] =
                    v->lattice_vectors[LATTICE_DIR_NW];
                index_voxel(ss, x - 1, y - 1)->lattice_vectors[LATTICE_DIR_SW] =
                    v->lattice_vectors[LATTICE_DIR_NE];

                // Keep track of stuff needed to plot force vector:
                #if 0
                barrierCount++;
                barrierxSum += x;
                barrierySum += y;
                barrierFx += nE[index] + nNE[index] + nSE[index] - nW[index] - nNW[index] - nSW[index];
                barrierFy += nN[index] + nNE[index] + nNW[index] - nS[index] - nSE[index] - nSW[index];
                #endif
            }
        }
    }
    #endif
}

void step_simulation_state(simulation_state_t* ss)
{
    // what happens if we set the boundary conditions every single step...?
    collide(ss);
    stream(ss);

    // TODO: move barriers

}


bool simulation_state_is_stable(simulation_state_t* ss)
{
    const int y = ss->params->ydim / 2;
    for (int x = 0; x < ss->params->xdim; x++) {
        fluid_voxel_t* v = index_voxel(ss, x, y);
        if (v->rho <= 0) return false;
    }

    return true;
}

int simulation_state_initialize_log_file(FILE* f, simulation_state_t* ss)
{
    fprintf(f, "curl: float[%5i, %5i]\r\n", ss->params->xdim, ss->params->ydim);
    fprintf(f, "density: float[%5i, %5i]\r\n", ss->params->xdim, ss->params->ydim);
    fprintf(f, "speed: float[%5i, %5i]\r\n", ss->params->xdim, ss->params->ydim);
    return 0;
}

static void compute_curl(simulation_state_t* ss)
{
    // interior sites only; leave edges set to zero
    for (int y = 1; y < (ss->params->ydim - 1); y++) {
        for (int x = 1; x < (ss->params->xdim - 1); x++) {
            fluid_voxel_t* v = index_voxel(ss, x, y);
            v->curl = (index_voxel(ss, x + 1, y)->u[1] -
                       index_voxel(ss, x - 1, y)->u[1] -
                       index_voxel(ss, x, y + 1)->u[0] +
                       index_voxel(ss, x, y - 1)->u[0]);
        }
    }
}

int simulation_state_append_sim_frame_to_log_file(FILE* f, simulation_state_t* ss)
{
    size_t dim = ss->params->xdim * ss->params->ydim;
    float* curl = malloc(dim * sizeof(float));
    float* density = malloc(dim * sizeof(float));
    float* speed = malloc(dim * sizeof(float));

    compute_curl(ss);

    int idx = 0;
    for (int y = 0; y < ss->params->ydim; y++) {
        for (int x = 0; x < ss->params->xdim; x++) {
            fluid_voxel_t* v = index_voxel(ss, x, y);
            curl[idx] = v->curl;
            density[idx] = v->rho;
            speed[idx] = sqrt((v->u[0] * v->u[0]) + (v->u[1] * v->u[1]));
            idx++;
        }
    }

    int succeeded = 1;
    succeeded &= (fwrite(curl, sizeof(float), dim, f) == dim);
    succeeded &= (fwrite(density, sizeof(float), dim, f) == dim);
    succeeded &= (fwrite(speed, sizeof(float), dim, f) == dim);

    return succeeded ? 0 : -1;
}
