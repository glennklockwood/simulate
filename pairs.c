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
           *restrict force = atoms->force, 
           *restrict stress = atoms->stress;
    double boxl2[MD_NDIMS];

    for ( int i = 0; i < MD_NDIMS; i++ )
        boxl2[i] = 2.0*system_info->boxli[i];

    neighbors = construct_neighbor_list( system_info, atoms );

    for ( int i = 0; i < atoms->count; i++ )
    {
        double pos_i[MD_NDIMS]; 
        for ( int k = 0; k < MD_NDIMS; k++ )
            pos_i[k] = pos[COORD(i,k)];

        /* 2-body pairs */
        for ( int m = neighbors->nlb[i]; m < neighbors->nle[i]; m++ )
        {
            double rij2 = 0.0;
            double potential_mag, force_mag;
            double rij[MD_NDIMS],
                   fij[MD_NDIMS];
            int j;
            j = neighbors->nlist[m];

            // calculate j-i interatomic distance
            for ( int k = 0; k < MD_NDIMS; k++ )
            {
                int mx;
                rij[k] = pos[COORD(j,k)] - pos_i[k];
                mx = (int)(rij[k] * boxl2[k]);
                rij[k] -= (double)mx * system_info->boxl[k];
            }

            for ( int k = 0; k < MD_NDIMS; k++ )
                rij2 += rij[k]*rij[k];

            if ( sqrt(rij2) > system_info->boxl[0] )
            {
                fprintf( stderr, "Interatomic spacing is larger than box!\n" );
                exit(1);
            }

            // calculate potential and force from this pair
            simulation_parameters->calculate_2body(
                atoms->type[i], 
                atoms->type[j],
                rij2, 
                simulation_parameters, 
                &potential_mag, 
                &force_mag);

            fprintf(stderr, "Interaction at %15g is %15g erg, %15g dyne\n",
                sqrt(rij2), potential_mag, force_mag * sqrt(rij2));

            for ( int k = 0; k < MD_NDIMS; k++ )
            {
                fij[k] = rij[k] * force_mag;
                /* still not sure these are going in the right direction */
                force[COORD(i,k)] += fij[k];
                force[COORD(j,k)] -= fij[k];
                stress[k] += fij[k] * rij[k];
/*              stress[k][k] += fij[k] * rij[k]; */
            }
        }
    }
//  dump_table( simulation_parameters );
    return;
}
