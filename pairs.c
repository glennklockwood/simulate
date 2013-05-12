/******************************************************************************
 *  @brief      Molecular Dynamics Simulation - Velocity Verlet Integrator
 *  @author     Glenn K. Lockwood
 *  @date       March 2013
 *
 *  @file       pairs.c
 ******************************************************************************/

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include "simulate.h"

/**
 * @brief analogous to pairs.f
 */
void calculate_new_forces
(  
    sim_params_t *simulation_parameters,
    system_info_t *system_info,
    atom_list_t *atoms 
)
{
    neighbor_list_verlet_t *neighbors;
    double *restrict pos = atoms->pos, 
           *restrict vel = atoms->vel, 
           *restrict accel = atoms->accel, 
           *restrict force = atoms->force, 
           *restrict stress = atoms->stress;
    double boxl2[MD_NDIMS];

    for ( int i = 0; i < MD_NDIMS; i++ )
        boxl2[i] = 2.0*system_info->boxli[i];

    neighbors = construct_neighbor_list( system_info, atoms );

    for ( int i = 0; i < atoms->count; i++ )
    {
        double pos_i[MD_NDIMS]; 
        for ( int m = 0; m < MD_NDIMS; m++ )
            pos_i[COORD(i,m)] = pos[COORD(i,m)];

        /* 2-body pairs */
        for ( int m = neighbors->nlb[i]; m < neighbors->nle[i]; m++ )
        {
            double rij2, potential_mag, force_mag;
            double rij[MD_NDIMS],
                   fij[MD_NDIMS];
            int j;
            j = neighbors->nlist[m];
            for ( int k = 0; k < MD_NDIMS; k++ )
            {
                rij[k] = pos[COORD(j,k)] - pos_i[k];
                rij[k] -= (double)((int)(rij[k] * boxl2[k])) * 
                    system_info->boxl[k];
            }

            simulation_parameters->calculate_2body(
                atoms->ltype[i], 
                atoms->ltype[j],
                rij2, 
                simulation_parameters, 
                &potential_mag, 
                &force_mag);

            for ( int k = 0; k < MD_NDIMS; k++ )
            {
                fij[k] = rij[k] * force_mag;
                force[COORD(i,k)] -= fij[k];
                force[COORD(j,k)] += fij[k];
                stress[k] += fij[k] * rij[k];
/*              stress[k][k] += fij[k] * rij[k]; */
            }
        }
    }
}

