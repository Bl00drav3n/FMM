/*
 * ============================================================================
 *
 *       Filename:  fmm_sim.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  07.02.2015 20:02:37
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Bl00drav3n (), 
 *   Organization:  
 *
 * ============================================================================
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "sim_math.h"
#include "memory.h"
#include "quad_tree.h"

struct fmm_data
{
    uint32_t num_particles;
    v2 *positions;
    v2 *velocities;
    v2 *forces;
    v2 dim;

    uint32_t num_levels;
    uint32_t num_coeffs;
    struct cell root_cell;
};

static struct fmm_data data = {};

static void
reset_forces()
{
    for(uint32_t i = 0; i < data.num_particles; i++) {
        data.forces[i] = V2(0,0);
    }
}

static inline v2
canonicalize(v2 v)
{
    v2 result;
    v2 half_dim = data.root_cell.half_dim;
    result.x = v.x - ((int32_t)(v.x / half_dim.x)) * data.dim.x;
    result.y = v.y - ((int32_t)(v.y / half_dim.y)) * data.dim.y;
    assert(result.x >= -half_dim.x &&
            result.x <= half_dim.x);
    assert(result.y >= -half_dim.y &&
            result.y <= half_dim.y);
    return result;
}

/* 
 * TODO: try to speed up the function (intrinsic?)
 * in general we don't want to call floorf
 */
static inline v2
canonicalize_position(v2 pos)
{
    v2 result;
    result.x = pos.x - 
        (floorf(pos.x / data.dim.x)) * data.dim.x;
    result.y = pos.y - 
        (floorf(pos.y / data.dim.y)) * data.dim.y;
    assert(result.x >= 0.f &&
            result.x < data.dim.x);
    assert(result.y >= 0.f &&
            result.y < data.dim.y);
    return result;
}

static inline v2
get_min_diff(v2 a, v2 b)
{
    v2 result;
    v2 diff = a - b;
    result = canonicalize(diff);
    return result;
}

static void
calculate_forces()
{
    /* O(n^2) approach */
    for(uint32_t j = 0; j < data.num_particles; j++) {
        v2 cur_pos = data.positions[j];
        for(uint32_t i = j + 1; i < data.num_particles; i++) {
            v2 other_pos = data.positions[i];
            v2 diff = get_min_diff(cur_pos, other_pos);
            float absForce = 1 / lengthSq(diff);
            v2 force = absForce * diff;
            data.forces[j] += force;
            data.forces[i] -= force;
        }
    }
}

static void
integrate(float timeStep)
{
    /* simple (inaccurate) integration */
    for(uint32_t i = 0; i < data.num_particles; i++) {
        v2 velocity = data.velocities[i];
        data.velocities[i] += timeStep * data.forces[i];
        data.positions[i] = 
            canonicalize_position(data.positions[i] + timeStep * velocity);
    }
}

static void
initialize_data(uint32_t memsize, 
        uint32_t num_particles_x, uint32_t num_particles_y, 
        uint32_t num_levels, uint32_t num_coeffs)
{
    uint32_t num_particles = num_particles_x * num_particles_y;
    memory.free_memory = (uint8_t*)calloc(memsize, 1);
    memory.free_memory_size = memsize;
    
    data.positions = push_array(v2, num_particles);
    data.velocities = push_array(v2, num_particles);
    data.forces = push_array(v2, num_particles);
    data.num_particles = num_particles;
    
    /* initialize equidistant grid */
    v2 min_corner = V2(0, 0);
    v2 max_corner = V2(1, 1);
    v2 dim = max_corner - min_corner;
    float dx = dim.x / num_particles_x;
    float dy = dim.y / num_particles_y;
    uint32_t idx = 0;
    for(uint32_t j = 0; j < num_particles_y; j++) {
        for(uint32_t i = 0; i < num_particles_x; i++) {
            data.positions[idx] =
                min_corner + V2(i * dx, j * dy);
            data.velocities[idx] = V2(0,0);
            idx++;
        }
    }

    data.dim = dim;
    data.num_levels = num_levels;
    data.num_coeffs = num_coeffs;

    v2 half_dim = 0.5f * data.dim;
    initialize_quad_tree(&data.root_cell, half_dim, half_dim, 
            num_levels, num_coeffs);
}

static void
debug_write_positions_to_file(const char *filename)
{
    FILE *file = fopen(filename, "w");
    if(file) {
        for(uint32_t i = 0; i < data.num_particles; i++) {
            fprintf(file, "%.3f %.3f\n", 
                    data.positions[i].x, data.positions[i].y);
        }
        fclose(file);
    } else {
        fprintf(stderr, "Error opening file %s\n", filename);
    }
}

/* 
 * NOTE: We call this on the root node, which will clear all the a of non-leaf
 * nodes and builds the multipole moments (a) for every remainging cell in 
 * the tree.
 */
static void
accumulate_cell_multipoles(struct cell *this_)
{
    if(this_ == 0)
        return;
    else if(this_->childs[0] != 0) {
        /* NOTE: clear a to (0,0) if this is not a leaf node. */
        for(uint32_t i = 0; i < data.num_coeffs; i++) {
            this_->a[i] = V2(0, 0);
        }
    }

    for(uint32_t i = 0; i < 4; i++) {
        accumulate_cell_multipoles(this_->childs[i]);
    }

    /* do the actual work */
    if(this_->parent != 0) {
        this_->parent->a[0] += this_->a[0];
        v2 z = this_->parent->center - this_->center;
        for(uint32_t k = 1; k < data.num_coeffs; k++) {
            this_->parent->a[k] -= (this_->a[0]/k) * pow(z, k);
            for(uint32_t l = 1; l <= k; l++) {
                this_->parent->a[k] += binomial(k - 1, l - 1) *
                    this_->a[l] * pow(z, k - l);
            }
        }
    }
}

static void
reset_leave_nodes(struct cell *this_)
{
    if(this_ == 0)
        return;
    else if(this_->childs[0] == 0) {
        for(uint32_t i = 0; i < data.num_coeffs; i++) {
            this_->a[i] = V2(0,0);
        }
    }

    for(uint32_t i = 0; i < 4; i++) {
        reset_leave_nodes(this_->childs[i]);
    }
}

/*
 * NOTE: When this function returns, we have multipole moments for all
 * cells in the tree.
 */
static void
calculate_multipoles()
{
    /* 
     * TODO: Clear all a to (0,0). We are doing this for all
     * non-leaf nodes in the upward pass already. We need a good way to do this
     * efficiently - maybe in the downward pass? For now we just use a little
     * helper function.
     */

    reset_leave_nodes(&data.root_cell);

    for(uint32_t i = 0; i < data.num_particles; i++) {
        v2 pos = data.positions[i];
        struct cell *residence = find_residence(&data.root_cell, pos);
        residence->a[0] += V2(1, 0);
        for(uint32_t k = 1; k < data.num_coeffs; k++) {
            /* TODO: make this faster (don't use pow, precompute 1/k) */
            residence->a[k] -= pow(pos - residence->center, k) / (float)k;
        }
    }

    accumulate_cell_multipoles(&data.root_cell);
}

int 
main(int argc, char *argv[])
{
    float timeStep = 0.1f;
    uint32_t nIntegrations = 1000;
    uint32_t num_particles_x = 16;
    uint32_t num_particles_y = 16;
    uint32_t num_particles = num_particles_x * num_particles_y;
    uint32_t num_levels = logf(num_particles) / logf(4.f);
    uint32_t num_coefficients = 4;
    uint32_t memory_size = 256 * 1024 * 1024;
    initialize_data(memory_size, num_particles_x, num_particles_y,
            num_levels, num_coefficients);

    for(uint32_t loop = 0; loop < nIntegrations; loop++) {
        calculate_multipoles();
        reset_forces();
        calculate_forces();
        integrate(timeStep);
    }

    const char filename[] = "out.txt";
    debug_write_positions_to_file(filename);

    return 0;
}

