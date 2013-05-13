/******************************************************************************
 *  @brief      Molecular Dynamics Simulation - Velocity Verlet Integrator
 *  @author     Glenn K. Lockwood
 *  @date       March 2013
 *
 *  @file       move.c
 ******************************************************************************/

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include "simulate.h"

/**
 * @brief velocity verlet integration algorithm
 */
void next_iteration
(
    sim_params_t *simulation_parameters,
    system_info_t *system_info,
    atom_list_t *atoms
)
{
    double *restrict pos = atoms->pos,
           *restrict vel = atoms->vel,
           *restrict accel = atoms->accel,
           *restrict force = atoms->force,
           *restrict inv_mass = simulation_parameters->time_over_mass,
           *restrict boxl = system_info->boxl,
           *restrict boxli = system_info->boxli;

    /* Velocity Verlet algorithm */
    for ( int i = 0; i < atoms->count; i++ )
    {
        double displacement = 0.0;
        for ( int j = 0; j < MD_NDIMS; j++ )
        {
            double t;
            int m, n;
            vel[COORD(i,j)] += accel[COORD(i,j)];
            pos[COORD(i,j)] += vel[COORD(i,j)];
            displacement += vel[COORD(i,j)] * vel[COORD(i,j)];

            t = pos[COORD(i,j)] * boxli[j];
            m = (int)t;
            n = (int)(t - 1.0);
            pos[COORD(i,j)] = (t - (double)(m + n)) * boxl[j];
        }
    }

    calculate_new_forces( simulation_parameters, system_info, atoms );

    for ( int i = 0; i < atoms->count; i++ )
    {
        for ( int j = 0; j < MD_NDIMS; j++ )
        {
            accel[COORD(i,j)] = force[COORD(i,j)] 
                * inv_mass[atoms->type[i]];
            vel[COORD(i,j)] += accel[COORD(i,j)];
        }
    }

    return;
}


