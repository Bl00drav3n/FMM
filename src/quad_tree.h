/*
 * =====================================================================================
 *
 *       Filename:  quad_tree.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08.02.2015 19:20:24
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Bl00drav3n (), 
 *   Organization:  
 *
 * =====================================================================================
 */

/* NOTE: Should we include a boolean is_leaf? (not neccessary) */
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

#include "quad_tree.cpp"
