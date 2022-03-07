/**
 * Code is adapted from the javascript found at https://physics.weber.edu/schroeder/fluids/
 * originally by Dan Schroeder
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "simulation.h"

static inline float clamp_f(float f, float min, float max)
{
    if (f < min) return min;
    if (f > max) return max;
    return f;
}

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
    v->lattice_vectors[LATTICE_DIR_E] = (1.0 / 9.0) * (1 + ux3       + 4.5*ux2        - u215);
    v->lattice_vectors[LATTICE_DIR_W] = (1.0 / 9.0) * (1 - ux3       + 4.5*ux2        - u215);
    v->lattice_vectors[LATTICE_DIR_N] = (1.0 / 9.0) * (1 + uy3       + 4.5*uy2        - u215);
    v->lattice_vectors[LATTICE_DIR_S] = (1.0 / 9.0) * (1 - uy3       + 4.5*uy2        - u215);
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
            printf("%6.4f ", (double)ss->voxels[(y * ss->params->xdim) + x].lattice_vectors[1]);
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

static void collide(simulation_state_t* restrict ss)
{
    const float omega = 1.0 / ((3 * ss->params->viscosity) + 0.5);      // reciprocal of relaxation time
    for (int y = 1; y < (ss->params->ydim - 1); y++) {
        for (int x = 1; x < (ss->params->xdim - 1); x++) {
            fluid_voxel_t* v = index_voxel(ss, x, y);

            //                             0   N   S   E   W  NE  SE  NW  SW
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

static void stream(simulation_state_t* restrict ss)
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
        for (int x = 1; x < (ss->params->xdim - 1); x++) {
            fluid_voxel_t* v = index_voxel(ss, x, y);

            v->lattice_vectors[LATTICE_DIR_W] = index_voxel(ss, x + 1, y)->lattice_vectors[LATTICE_DIR_W];
            v->lattice_vectors[LATTICE_DIR_SW] = index_voxel(ss, x + 1, y + 1)->lattice_vectors[LATTICE_DIR_SW];
        }
    }

    // Now handle bounce-back from barriers
    for (int b = 0; b < ss->num_barriers; b++) {
        barrier_t* barrier = &ss->barriers[b];

        float F_barrier[2] = { 0 };

        for (int by = 0; by < barrier->ydim; by++) {
            for (int bx = 0; bx < barrier->xdim; bx++) {
                if (barrier->occupancy[by * barrier->xdim + bx]) {
                    int x = bx + (int)barrier->x;
                    int y = by + (int)barrier->y;

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
                    F_barrier[0] += (v->lattice_vectors[LATTICE_DIR_E] +
                                     v->lattice_vectors[LATTICE_DIR_NE] +
                                     v->lattice_vectors[LATTICE_DIR_SE] -
                                     v->lattice_vectors[LATTICE_DIR_W] -
                                     v->lattice_vectors[LATTICE_DIR_NW] -
                                     v->lattice_vectors[LATTICE_DIR_SW]);

                    F_barrier[1] += (v->lattice_vectors[LATTICE_DIR_N] +
                                     v->lattice_vectors[LATTICE_DIR_NE] +
                                     v->lattice_vectors[LATTICE_DIR_NW] -
                                     v->lattice_vectors[LATTICE_DIR_S] -
                                     v->lattice_vectors[LATTICE_DIR_SE] -
                                     v->lattice_vectors[LATTICE_DIR_SW]);
                }
            }
        }

        // update the barrier's velocity
        float k_force[2] = {
            -(barrier->x - barrier->anchor[0]) * barrier->k[0],
            -(barrier->y - barrier->anchor[1]) * barrier->k[1]
        };

        printf("f_barrier = {%9.2f, %9.2f}\n", F_barrier[0], F_barrier[1]);
        printf("f_spring  = {%9.2f, %9.2f}\n", k_force[0], k_force[1]);
        printf("barrier_pos = {%9.2f, %9.2f}\n\n", barrier->x, barrier->y);

        if (barrier->k[0] < 1e6)
            barrier->u[0] += (F_barrier[0] + k_force[0]) / barrier->mass;
        if (barrier->k[1] < 1e6)
            barrier->u[1] += (F_barrier[1] + k_force[1]) / barrier->mass;
    }
}

static void move_barriers(simulation_state_t* restrict ss)
{
    for (int b = 0; b < ss->num_barriers; b++) {
        barrier_t* barrier = &ss->barriers[b];

        barrier->x = clamp_f(barrier->x + barrier->u[0], 0, ss->params->xdim - barrier->xdim);
        barrier->y = clamp_f(barrier->y + barrier->u[1], 0, ss->params->ydim - barrier->ydim);

        // TODO: push fluid
    }
}

void step_simulation_state(simulation_state_t* restrict ss)
{
    move_barriers(ss);

    // what happens if we set the boundary conditions every single step...?
    collide(ss);
    stream(ss);

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

/**
 * File format:
 * A header which specifies the fields in the raw data section.
 * Each line of the header should be terminated by a semicolon.
 * The header ends with the special string "FSIM".
 *
 * All data after the letters "FSIM" are packed binary frames that consist of the
 *    header:
 *        "<field_0_name>: <field_0_datatype>[<field_0_dimensions>];"
 *        "<field_1_name>: <field_1_datatype>[<field_1_dimensions>];"
 *        ....
 *        "<field_n-1_name>: <field_n-1_datatype>[<field_n-1_dimensions>];"
 *
 *    data ():
 *        "FSIM" frame_0_field_0_data frame_0_field_1_data ... frame_0_field_n-1_data
 *               frame_1_field_0_data frame_1_field_1
 */
int simulation_state_initialize_log_file(FILE* f, simulation_state_t* ss)
{
    fprintf(f, "curl: f[%5i, %5i];\r\n", ss->params->xdim, ss->params->ydim);
    fprintf(f, "density: f[%5i, %5i];\r\n", ss->params->xdim, ss->params->ydim);
    fprintf(f, "speed: f[%5i, %5i];\r\n", ss->params->xdim, ss->params->ydim);
    fprintf(f, "barriers: ?[%5i, %5i];\r\n", ss->params->xdim, ss->params->ydim);
    fprintf(f, "FSIM");

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

    uint8_t* barriers = calloc(dim, sizeof(uint8_t));
    for (int i = 0; i < ss->num_barriers; i++) {
        const barrier_t* barr = &ss->barriers[i];
        for (int by = 0; by < barr->ydim; by++) {
            for (int bx = 0; bx < barr->xdim; bx++) {
                int x = barr->x + bx;
                int y = barr->y + by;

                barriers[(y * ss->params->xdim) + x] |= barr->occupancy[by * barr->xdim + bx];
            }
        }
    }

    int succeeded = 1;
    succeeded &= (fwrite(curl, sizeof(float), dim, f) == dim);
    succeeded &= (fwrite(density, sizeof(float), dim, f) == dim);
    succeeded &= (fwrite(speed, sizeof(float), dim, f) == dim);
    succeeded &= (fwrite(barriers, sizeof(uint8_t), dim, f) == dim);

    return succeeded ? 0 : -1;
}

barrier_t* barrier_create_manual(int xdim, int ydim, int x, int y, const bool* occupancy)
{
    barrier_t* barr = calloc(1, sizeof(barrier_t));

    barr->xdim = xdim;
    barr->ydim = ydim;

    barr->x = x;
    barr->y = y;
    barr->anchor[0] = x;
    barr->anchor[1] = y;
    barr->k[0] = 1e8;
    barr->k[1] = .1;
    barr->mass = 1000;
    barr->occupancy = calloc(xdim * ydim, sizeof(bool));
    memcpy(barr->occupancy, occupancy, xdim * ydim * sizeof(bool));

    return barr;
}

barrier_t* barrier_create_rectangle(int width, int height)
{
    barrier_t* barr = NULL;
    return barr;
}

void barrier_destroy(barrier_t* barr)
{
    free(barr->occupancy);
    free(barr);
}

void simulation_state_add_barrier(simulation_state_t* ss, const barrier_t* barr)
{
    ss->num_barriers++;
    ss->barriers = realloc(ss->barriers, ss->num_barriers * sizeof(barrier_t));

    memcpy(&ss->barriers[ss->num_barriers - 1], barr, sizeof(barrier_t));
}
