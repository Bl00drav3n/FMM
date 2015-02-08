/*
 * =====================================================================================
 *
 *       Filename:  memory.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08.02.2015 19:21:08
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Bl00drav3n (), 
 *   Organization:  
 *
 * =====================================================================================
 */

struct memory_arena
{
    uint8_t *free_memory;
    uint32_t free_memory_size;
};

#include "memory.cpp"

#define push_struct(type) (type*)push_memory(sizeof(type))
#define push_array(type, amount) (type*)push_memory(amount * sizeof(type))
