/******************************************************************************
 *  @brief      Molecular Dynamics Simulation
 *  @author     Glenn K. Lockwood
 *  @date       March 2013
 *
 *  @file       simulate.c
 ******************************************************************************/

/** @mainpage   simulate
 *              This is an experimental molecular dynamics simulation code 
 *              developed to resemble the MOLDYN code developed at the 
 *              Interfacial Molecular Science Laboratory at Rutgers University.
 *              Rewritten in C99, it strives to function similarly to the 
 *              original FORTRAN77 code while allowing for the integration of
 *              newer, high-performance libraries and APIs such as PAPI, 
 *              OpenACC, and stronger automatic vectorization support.
 *              
 ******************************************************************************/

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include "simulate.h"

/******************************************************************************
 *
 *   Functions
 *
 ******************************************************************************/
int main( int argc, char **argv )
{
    sim_params_t *simulation_parameters;
    system_info_t *system_info;
    atom_list_t *atoms;
    
    simulation_parameters = malloc(sizeof(*simulation_parameters));
    system_info = malloc(sizeof(*system_info));
    atoms = malloc(sizeof(*atoms));
    if ( !simulation_parameters || !system_info || !atoms )
        malloc_error ();

    init_simulation( simulation_parameters, system_info );
    get_input_params( simulation_parameters );
    atoms = get_starting_config( simulation_parameters, system_info );
    modify_starting_config( system_info, atoms );

    simulate( simulation_parameters, system_info, atoms );

    return 0;
}

/** 
 * @brief perform the initial memory allocation for simulation_parameters and 
 * system_info.  requires knowledge of the number of elements 
 */
void init_simulation
(
    sim_params_t *simulation_parameters,
    system_info_t *system_info 
)
{
    if ( !(simulation_parameters = malloc(sizeof(*simulation_parameters)))
    ||   !(system_info = malloc(sizeof(*system_info)) ) )
        malloc_error();

    return;
}

/**
 * @brief function to control the calculation the trajectory of atoms
 */
void simulate
(
    sim_params_t *simulation_parameters,
    system_info_t *system_info,
    atom_list_t *atoms 
)
{
    bool complete = false;
    while ( !complete )
    {
        next_iteration( simulation_parameters, system_info, atoms );
    }
}

/**
 * @brief analogous to move.f
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
           *restrict stress = atoms->stress,
           *restrict inv_mass = simulation_parameters->inv_mass,
           *restrict boxl = system_info->boxl,
           *restrict boxli = system_info->boxli;

    /* Velocity Verlet algorithm */
    for ( int i = 0; i < atoms->count; i++ )
    {
        for ( int j = 0; j < MD_NDIMS; j++ )
        {
            double t;
            int m, n;
            vel[COORD(i,j)] += accel[COORD(i,j)];
            pos[COORD(i,j)] += vel[COORD(i,j)];

            t = pos[COORD(i,j)] * boxli[j];
            m = (int)t;
            n = m - 1;
            pos[COORD(i,j)] = (t - (double)(m - n)) * boxl[j];
        }
    }

    calculate_new_forces( simulation_parameters, system_info, atoms );

    for ( int i = 0; i < atoms->count; i++ )
    {
        for ( int j = 0; j < MD_NDIMS; j++ )
        {
            accel[COORD(i,j)] = force[COORD(i,j)] * inv_mass[i];
            vel[COORD(i,j)] += accel[COORD(i,j)];
        }
    }

    return;
}

/**
 * @brief throw the appropriate error for a memory allocation failure
 */
void malloc_error( void )
{
    fprintf(stderr, "FATAL ERROR: could not allocate memory\n");
    exit(MD_ERR_MALLOC);
}

/**
 * @brief do any pre-simulation system manipulation.  from aamoldyn.f
 */
void modify_starting_config( system_info_t *system_info, atom_list_t *atoms )
{
    return;
}

/**
 * @brief read in simulation parameters from IMSL tape5 file.  from aamoldyn.f
 */
void get_input_params( sim_params_t *simulation_parameters )
{
    return;
}
 
/**
 * @brief load an initial configuration in IMSL's knite12 format.  from aamoldyn.f
 */
atom_list_t *get_starting_config
(
    sim_params_t *simulation_parameters,
    system_info_t *system_info 
)
{
    FILE *fp;
    int garbage_int; 
    size_t size_int, size_dble, size_dble_ptr, size_array;
    int i, j, nbin;
    int reclength;
    int bytesread = 0;
    double garbage_dble;
    atom_list_t *atoms;

    size_int = sizeof(int);
    size_dble = sizeof(double);
    size_dble_ptr = sizeof(double*);

    if ( (fp = fopen( "knite12", "r" )) == NULL )
    {
        fprintf(stderr, "FATAL ERROR: could not open input file (knite12)\n");
        exit(MD_ERR_FOPEN);
    }

    if ( (atoms = malloc(sizeof(*atoms))) == NULL )
        malloc_error();
    
    rewind(fp);

    fread( &reclength, sizeof(reclength), 1, fp );

// read in the header with the pertinent dimension info 
    fread( &garbage_dble, size_dble, 1, fp); 
    fread( &(atoms->count), size_int, 1, fp);
    fread( &garbage_int, size_int, 1, fp);
    fread( &garbage_int, size_int, 1, fp);
    fread( &(simulation_parameters->num_ltypes), size_int, 1, fp);
    fread( &garbage_int, size_int, 1, fp);
    fread( &nbin, size_int, 1, fp);
    bytesread += size_dble + size_int*6;

// allocate all the necessary memory
    size_array = MD_NDIMS * atoms->count * size_dble;
    if ( (atoms->pos      = malloc(size_array)) == NULL
    ||   (atoms->vel      = malloc(size_array)) == NULL
    ||   (atoms->accel    = malloc(size_array)) == NULL
    ||   (atoms->force    = malloc(size_array)) == NULL
    ||   (atoms->ltype    = malloc(atoms->count * size_int)) == NULL
    ||   (atoms->index    = malloc(atoms->count * size_int)) == NULL )
        malloc_error();

    for ( i = 0; i < atoms->count; i++ )
        atoms->index[i] = i;

    // old BMH parameters stored in knite9
    for ( i = 0; i < simulation_parameters->num_ltypes; i++ )
    {
        for ( j = 0; j < simulation_parameters->num_ltypes; j++ )
        {
            fread(&garbage_dble,    size_dble, 1, fp );
            fread(&garbage_dble,    size_dble, 1, fp );
            fread(&garbage_dble,    size_dble, 1, fp );
            fread(&garbage_dble,    size_dble, 1, fp );
            fread(&garbage_dble,    size_dble, 1, fp );
            fread(&garbage_dble,    size_dble, 1, fp );
        }
        fread(&garbage_dble,size_dble,1,fp);
        fread(&garbage_dble,size_dble,1,fp);
    }
    bytesread += (6 * simulation_parameters->num_ltypes+ 2) 
                * simulation_parameters->num_ltypes * size_dble;

    // read in positions, velocities, accelerations
    // 18 = positions and 5 derivatives (6 total) in 3 dimensions (6x3=18)
    for ( i = 0; i < atoms->count; i++ )
    {
        double atom_vector[18];
        fread( atom_vector, size_dble, 18, fp );
        for ( j = 0; j < MD_NDIMS; j++ )
        {
            atoms->pos[COORD(i,j)] = atom_vector[j];
            atoms->vel[COORD(i,j)] = atom_vector[MD_NDIMS+j];
            atoms->accel[COORD(i,j)] = atom_vector[2*MD_NDIMS+j];
            atoms->force[COORD(i,j)] = 0.0;
        }
    }
    bytesread += atoms->count * 18 * size_dble;

    for ( j = 0; j < simulation_parameters->num_ltypes; j++ )
    {
      for ( i = 0; i < nbin; i++ )
      {
        fread(&garbage_dble,size_dble, 1, fp );
        fread(&garbage_int,size_int, 1, fp );
      }
    }
    bytesread += simulation_parameters->num_ltypes * nbin * (size_dble+size_int);

    fread( &(system_info->boxl), size_dble, 3, fp );
    fread( &garbage_dble, size_dble, 1, fp );
    fread( &garbage_dble, size_dble, 1, fp );
    fread( &garbage_dble, size_dble, 1, fp );
    bytesread += 6*size_dble;

    for ( i = 0; i < atoms->count; i++ )
    {
      fread( &(atoms->ltype[i]), size_int, 1, fp );
      fread( &garbage_int, size_int, 1, fp );
      fread( &garbage_int, size_int, 1, fp );
      fread( &garbage_int, size_int, 1, fp );
      atoms->ltype[i]--;
    }
    bytesread += 4 * atoms->count * size_int;

    fread( &reclength, sizeof(reclength), 1, fp );
    if ( reclength != bytesread )
    {
        fprintf( stderr, "ERROR: knite9 read in %d bytes of %d[%x] expected\n",
            bytesread, reclength, reclength );
        exit( MD_ERR_INPUT_FMT );
    }
    return atoms;
}
