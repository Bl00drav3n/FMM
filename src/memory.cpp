/*
 * =====================================================================================
 *
 *       Filename:  memory.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08.02.2015 19:22:14
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Bl00drav3n (), 
 *   Organization:  
 *
 * =====================================================================================
 */

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

