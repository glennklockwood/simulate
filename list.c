/******************************************************************************
 *  @brief      Molecular Dynamics Simulation - Verlet Neighbor List Generator
 *  @author     Glenn K. Lockwood
 *  @date       March 2013
 *
 *  @file       list.c
 ******************************************************************************/

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include "simulate.h"

/**
 * @brief analogous to list.f
 *
 * @note This subroutine is not implemented yet.
 *
 */
neighbor_list_verlet_t *construct_neighbor_list
(  
    system_info_t *system_info,
    atom_list_t *atoms 
)
{
    static neighbor_list_verlet_t *neighbor_list;

    create_link_cell( system_info, atoms );

    if ( neighbor_list )
    {
        free(neighbor_list->nlist);
        free(neighbor_list->nlb);
        free(neighbor_list->nle);
        free(neighbor_list);
    }

//  fprintf( stderr, "neighbor list with %d atoms\n", atoms->count );
    neighbor_list = malloc(sizeof(*neighbor_list));
    neighbor_list->nlist = malloc(sizeof(int)*atoms->count*MAX_NEIGHBORS);
    neighbor_list->nlb = malloc(sizeof(int)*atoms->count);
    neighbor_list->nle = malloc(sizeof(int)*atoms->count);
    if ( !neighbor_list 
    ||   !neighbor_list->nlist 
    ||   !neighbor_list->nlb 
    ||   !neighbor_list->nle )
        malloc_error();

    /*
     * populate the neighbor list here
     */
    int frame = 0;
    for ( int i = 0; i < atoms->count-1; i++ )
    {
        neighbor_list->nlb[i] = frame;
        for ( int j = i+1; j < atoms->count; j++ )
            neighbor_list->nlist[frame++] = j;
        neighbor_list->nle[i] = frame;
    }
    neighbor_list->nlb[atoms->count-1] = frame;
    neighbor_list->nle[atoms->count-1] = frame;

    return neighbor_list;
}
