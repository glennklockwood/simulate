/******************************************************************************
 *  @brief      Molecular Dynamics Simulation - Potential Energy Functions
 *  @author     Glenn K. Lockwood
 *  @date       March 2013
 *
 *  @file       table.c
 ******************************************************************************/

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include "simulate.h"

#define LJ_EPSILON 0
#define LJ_SIGMA 1
/**
 * @brief analogous to the routines in table.f
 */
void pair_lj
( 
    int itype,
    int jtype,
    double rij2, 
    sim_params_t *simulation_parameters,
    double *potential_mag, 
    double *force_mag 
)
{
    double rij6 = rij2*rij2*rij2;
    double rij12 = rij6*rij6;

    double epsilon4, sigma, sigma12, sigma6;
    double part1, part2;

    int kt = LT2KT(itype,jtype,simulation_parameters->num_ltypes);

    epsilon4 = 4.0 * simulation_parameters->pair_param[kt][LJ_EPSILON];
    sigma = simulation_parameters->pair_param[kt][LJ_SIGMA];
    sigma6 = sigma*sigma;
    sigma6 *= sigma6*sigma6;
    sigma12 = sigma6*sigma6;

    part1 = epsilon4 * sigma12 / rij12;
    part2 = epsilon4 * sigma6 / rij6;

    *potential_mag = part1 - part2;
    *force_mag = part2 - 2.0*part1;
    *force_mag *= 6.0 / sqrt(rij2);

    return;
}
#undef LJ_EPSILON
#undef LJ_SIGMA


