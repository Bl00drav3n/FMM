/*
 * =====================================================================================
 *
 *       Filename:  quad_tree.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08.02.2015 19:25:47
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Bl00drav3n (), 
 *   Organization:  
 *
 * =====================================================================================
 */

static void
initialize_cell(struct cell *this_, struct cell *parent, 
        v2 center, v2 half_dim, uint32_t num_coeffs)
{
    this_->parent = parent;
    this_->half_dim = half_dim;
    this_->center = center;
    this_->a = push_array(v2, num_coeffs);
    this_->b = push_array(v2, num_coeffs);
}

static void
create_children(struct cell *this_, uint32_t num_coeffs)
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
        initialize_cell(cur_child, this_, centers[i], hhalf_dim, num_coeffs);
        this_->childs[i] = cur_child;
    }
}

static void
subdivide(struct cell *this_, uint32_t level, uint32_t num_coeffs)
{
    if(level == 1) {
        return;
    }

    create_children(this_, num_coeffs);
    for(uint32_t i = 0; i < 4; i++) {
        subdivide(this_->childs[i], level - 1, num_coeffs);
    }
}

static void
initialize_quad_tree(struct cell *root, 
        v2 center, v2 half_dim, uint32_t num_levels, uint32_t num_coeffs)
{
    initialize_cell(root, 0, center, half_dim, num_coeffs);
    subdivide(root, num_levels, num_coeffs);
}

