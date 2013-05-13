/******************************************************************************
 *  @brief      Molecular Dynamics Simulation - Link-Cell Method
 *  @author     Glenn K. Lockwood
 *  @date       March 2013
 *
 *  @file       cells.c
 ******************************************************************************/

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include "simulate.h"

/**
 * @brief analogous to cells.f
 *
 * @note This subroutine is incomplete and needs to be finished before this
 *       code has a chance at working.
 */
cell_map_t *create_link_cell( system_info_t *system_info, atom_list_t *atoms )
{
    static cell_map_t *cell_list; 

    if ( cell_list )
    {
        free(cell_list);
    }

    cell_list = malloc(sizeof(*cell_list));

    if (!cell_list)
        malloc_error();

    /*
     *  maxcel = maxn
     *  mapsiz = 35*maxcel
     *
     *  lhead(0:maxcel)
     */

   /*
    cell_list->
    cell_list->lhead = malloc(sizeof(int)*cell_list->ncell);
    cell_list->llist = malloc(sizeof(int)*cell_list->ncell);
    cell_list->map
    */

    return cell_list;
}


