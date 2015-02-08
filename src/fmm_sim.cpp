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

struct memory_arena
{
    uint8_t *free_memory;
    uint32_t free_memory_size;
};

static struct memory_arena memory = {};

static void *
push_memory(uint32_t size)
{
    void *result;
    assert(memory.free_memory_size >= size);
    result = memory.free_memory;
    memory.free_memory += size;
    memory.free_memory_size -= size;
    return result;
}

#define push_struct(type) (type*)push_memory(sizeof(type))
#define push_array(type, amount) (type*)push_memory(amount * sizeof(type))

struct cell
{
    v2 center;
    v2 half_dim;
    struct cell *parent;
    union {
        struct cell *childs[4];
        struct cell *childs2d[2][2];
    };

    v2 *a;
    v2 *b;
};

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
initialize_cell(struct cell *this_, struct cell *parent, 
        v2 center, v2 half_dim)
{
    this_->parent = parent;
    this_->half_dim = half_dim;
    this_->center = center;
    this_->a = push_array(v2, data.num_coeffs);
    this_->b = push_array(v2, data.num_coeffs);
}

static void
create_children(struct cell *this_)
{
    v2 centers[4];
    v2 center = this_->center;
    v2 hhalf_dim = 0.5f * this_->half_dim;
    centers[0] = center + V2(-hhalf_dim.x, -hhalf_dim.y);
    centers[1] = center + V2(hhalf_dim.x, -hhalf_dim.y);
    centers[2] = center + V2(-hhalf_dim.x, hhalf_dim.y);
    centers[3] = center + V2(hhalf_dim.x, hhalf_dim.y);

    struct cell *childs = push_array(struct cell, 4);
    for(uint32_t i = 0; i < 4; i++) {
        struct cell *cur_child = childs + i;
        initialize_cell(cur_child, this_, centers[i], hhalf_dim);
        this_->childs[i] = cur_child;
    }
}


static void
subdivide(struct cell *this_, uint32_t level)
{
    if(level == 1) {
        return;
    }

    create_children(this_);
    for(uint32_t i = 0; i < 4; i++) {
        subdivide(this_->childs[i], level - 1);
    }
}

static void
initialize_quad_tree(struct cell *root)
{
    v2 half_dim = 0.5f * data.dim;
    initialize_cell(root, 0, half_dim, half_dim);
    /* recursion */
    subdivide(root, data.num_levels);
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

    initialize_quad_tree(&data.root_cell);
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

int 
main(int argc, char *argv[])
{
    float timeStep = 0.1f;
    uint32_t nIntegrations = 1000;
    uint32_t num_particles_x = 16;
    uint32_t num_particles_y = 16;
    uint32_t num_particles = num_particles_x * num_particles_y;
    uint32_t num_levels = 4;
    uint32_t num_coefficients = 4;
    uint32_t memory_size = 256 * 1024 * 1024;
    initialize_data(memory_size, num_particles_x, num_particles_y,
            num_levels, num_coefficients);

    for(uint32_t loop = 0; loop < nIntegrations; loop++) {
        reset_forces();
        calculate_forces();
        integrate(timeStep);
    }

    const char filename[] = "out.txt";
    debug_write_positions_to_file(filename);

    return 0;
}

