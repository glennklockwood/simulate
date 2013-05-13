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
    
    init_simulation( &simulation_parameters, &system_info );
    get_input_params( simulation_parameters );
/*  atoms = get_starting_config( simulation_parameters, system_info ); */
    atoms = get_starting_config_xyz( simulation_parameters, system_info );
    init_constants( simulation_parameters, system_info );
    modify_starting_config( system_info, atoms );

    simulate( simulation_parameters, system_info, atoms );

    return 0;
}

/** 
 * @brief perform the initial memory allocation for simulation_parameters and 
 * system_info.  requires knowledge of the number of ltypes 
 */
void init_simulation
(
    sim_params_t **simulation_parameters,
    system_info_t **system_info 
)
{
    if ( !(*simulation_parameters = malloc(sizeof(**simulation_parameters)))
    ||   !(*system_info = malloc(sizeof(**system_info)) ) )
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
    int iterations = 0;
    bool complete = false;
    while ( !complete )
    {
        next_iteration( simulation_parameters, system_info, atoms );
        if ( ++iterations % 100 == 0 )
        {
            FILE *fp;
            int i, j;
            if ( (fp = fopen( "output.xyz", "a" )) )
            {
                fprintf( fp, "%d\n\n", atoms->count );
                for ( i = 0; i < atoms->count; i++ )
                {
                    int type = atoms->type[i];
                    fprintf(fp, "%d ", simulation_parameters->ltype[type]);
                    for ( j = 0; j < MD_NDIMS; j++ )
                        fprintf( fp, "%g ", atoms->pos[COORD(i,j)]*1e8 );
                    fprintf(fp, "\n");
                }
                fclose(fp);
                fprintf( stderr, "Wrote output at step %d\n", iterations );
            }
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
void get_input_params( sim_params_t *sp )
{
    int nlt, npair, nparam;
    /* 
     * currently a lot of hard-coded garbage for testing 
     */
    sp->num_types = 1;
    nlt = sp->num_types;
    npair = LT2KT(nlt-1, nlt-1, nlt) + 1;
    nparam = 2;
    fprintf( stderr, "Num types=%d\nNum pairs=%d\n", nlt, npair );

    /* number of atom types */
    sp->ltype = malloc(sizeof(*(sp->ltype)) * nlt);
    sp->mass = malloc(sizeof(*(sp->mass)) * nlt);
    sp->time_over_mass = malloc(sizeof(*(sp->time_over_mass)) * nlt);
    if ( !sp->ltype ) malloc_error();
    if ( !sp->mass ) malloc_error();
    if ( !sp->time_over_mass ) malloc_error();

    /* initialize everything */
    for ( int i = 0; i < nlt; i++ )
    {
        sp->ltype[i] = -1;
        sp->mass[i] = 0.0;
        sp->time_over_mass[i] = 0.0;
    }

    /* number of pair types */
    sp->pair_param = malloc(sizeof(*sp->pair_param) * npair);
    if ( !sp->pair_param ) malloc_error();
    for ( int i = 0; i < npair; i++ )
    {
        sp->pair_param[i] = malloc(sizeof(*sp->pair_param[i])*nparam);
        if ( !sp->pair_param[i] ) malloc_error();
        /* initialize everything to zero */
        for ( int j = 0; j < nparam; j++ )
            sp->pair_param[i][j] = 0.0;
    }

    /*
     * hard-coded garbage
     */
    sp->ltype[0] = 1;
    sp->mass[0] = 39.948 / 6.022e23;    // argon
    sp->timestep = 1.0e-15;
    sp->calculate_2body = &pair_lj;
    sp->calculate_3body = NULL;
    sp->pair_param[0][0] = 1.65e-14;    // 1.65e-21 J in ergs
    sp->pair_param[0][1] = 3.4e-08;     // cm

    return;
}
 
/**
 * @brief load an initial configuration in IMSL's knite12 format.  from 
 * aamoldyn.f
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
    fread( &(simulation_parameters->num_types), size_int, 1, fp);
    fread( &garbage_int, size_int, 1, fp);
    fread( &nbin, size_int, 1, fp);
    bytesread += size_dble + size_int*6;

// allocate all the necessary memory
    size_array = MD_NDIMS * atoms->count * size_dble;
    if ( (atoms->pos      = malloc(size_array)) == NULL
    ||   (atoms->vel      = malloc(size_array)) == NULL
    ||   (atoms->accel    = malloc(size_array)) == NULL
    ||   (atoms->force    = malloc(size_array)) == NULL
    ||   (atoms->stress   = malloc(size_array)) == NULL
    ||   (atoms->type     = malloc(atoms->count * size_int)) == NULL
    ||   (atoms->index    = malloc(atoms->count * size_int)) == NULL )
        malloc_error();

    for ( i = 0; i < atoms->count; i++ )
        atoms->index[i] = i;

    // old BMH parameters stored in knite9
    for ( i = 0; i < simulation_parameters->num_types; i++ )
    {
        for ( j = 0; j < simulation_parameters->num_types; j++ )
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
    bytesread += (6 * simulation_parameters->num_types+ 2) 
                * simulation_parameters->num_types * size_dble;

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
            atoms->stress[COORD(i,j)] = 0.0;
        }
    }
    bytesread += atoms->count * 18 * size_dble;

    for ( j = 0; j < simulation_parameters->num_types; j++ )
    {
      for ( i = 0; i < nbin; i++ )
      {
        fread(&garbage_dble,size_dble, 1, fp );
        fread(&garbage_int,size_int, 1, fp );
      }
    }
    bytesread += simulation_parameters->num_types * nbin * (size_dble+size_int);

    fread( &(system_info->boxl), size_dble, 3, fp );
    fread( &garbage_dble, size_dble, 1, fp );
    fread( &garbage_dble, size_dble, 1, fp );
    fread( &garbage_dble, size_dble, 1, fp );
    bytesread += 6*size_dble;

    for ( i = 0; i < atoms->count; i++ )
    {
      fread( &(atoms->type[i]), size_int, 1, fp );
      fread( &garbage_int, size_int, 1, fp );
      fread( &garbage_int, size_int, 1, fp );
      fread( &garbage_int, size_int, 1, fp );
      atoms->type[i]--;
    }
    bytesread += 4 * atoms->count * size_int;

    fread( &reclength, sizeof(reclength), 1, fp );
    if ( reclength != bytesread )
    {
        fprintf( stderr, "ERROR: knite9 read in %d bytes of %d[%x] expected\n",
            bytesread, reclength, reclength );
        exit( MD_ERR_INPUT_FMT );
    }
    fclose(fp);

    return atoms;
}

/**
 * @brief load an initial configuration in xyz format
 */
atom_list_t *get_starting_config_xyz
(
    sim_params_t *simulation_parameters,
    system_info_t *system_info 
)
{
    FILE *fp;
    size_t size_int, size_dble, size_dble_ptr, size_array;
    int i, j;
    atom_list_t *atoms;

    size_int = sizeof(int);
    size_dble = sizeof(double);
    size_dble_ptr = sizeof(double*);

    if ( (fp = fopen( "knite12.xyz", "r" )) == NULL )
    {
        fprintf(stderr, "FATAL ERROR: could not open input file (knite12)\n");
        exit(MD_ERR_FOPEN);
    }

    if ( (atoms = malloc(sizeof(*atoms))) == NULL )
        malloc_error();
    
    rewind(fp);

    fscanf(fp, "%d", &(atoms->count));
    fprintf( stderr, "Expecting %d atoms total\n", atoms->count );

// allocate all the necessary memory
    size_array = MD_NDIMS * atoms->count * size_dble;
    if ( (atoms->pos      = malloc(size_array)) == NULL
    ||   (atoms->vel      = malloc(size_array)) == NULL
    ||   (atoms->accel    = malloc(size_array)) == NULL
    ||   (atoms->force    = malloc(size_array)) == NULL
    ||   (atoms->stress   = malloc(size_array)) == NULL
    ||   (atoms->type     = malloc(atoms->count * size_int)) == NULL
    ||   (atoms->index    = malloc(atoms->count * size_int)) == NULL )
        malloc_error();

    for ( i = 0; i < atoms->count; i++ )
        atoms->index[i] = i;

    // read in positions
    for ( i = 0; i < atoms->count; i++ )
    {
        fscanf(fp, "%d", &(atoms->type[i]) );
        atoms->type[i]--;
/*** 
 *** TODO: fully separate ltype from type; have type be determined 
 ***       automatically by input file format 
 ***/
        for ( j = 0; j < MD_NDIMS; j++ )
        {
            fscanf(fp,"%lf", &(atoms->pos[COORD(i,j)]));
            atoms->vel[COORD(i,j)] = 0.0;
            atoms->accel[COORD(i,j)] = 0.0;
            atoms->force[COORD(i,j)] = 0.0;
        }
        fprintf( stderr, "Got atom of type %2d and pos ", atoms->type[i] );
        for ( j = 0; j < MD_NDIMS; j++ )
            fprintf( stderr, "%g ", atoms->pos[COORD(i,j)] );
        fprintf( stderr, "\n");

    }

    fclose(fp);

    // estimate the box length
    double maxl[MD_NDIMS];
    for ( j = 0; j < MD_NDIMS; j++ )
        maxl[j] = 0.0;

    for ( i = 0; i < atoms->count; i++ )
        for ( j = 0; j < MD_NDIMS; j++ )
            if ( atoms->pos[COORD(i,j)] > maxl[j] )
                maxl[j] = atoms->pos[COORD(i,j)];

    for ( j = 0; j < MD_NDIMS; j++ )
    {
        system_info->boxl[j] = maxl[j];
        system_info->boxl[j] *= 1.1;
    }

    fprintf( stderr, "Got %d atoms\nBox dims: ",
        atoms->count);
    for ( j = 0; j < MD_NDIMS; j++ )
        fprintf( stderr, "%g ", system_info->boxl[j] );
    fprintf( stderr, "\n" );

    return atoms;
}

/**
 * @brief initialize internal constants based on inputted data
 */
void init_constants
(
    sim_params_t *simulation_parameters,
    system_info_t *system_info 
)
{
    // inverse box dimensions
    for ( int i = 0; i < MD_NDIMS; i++ )
    {
        if ( system_info->boxl[i] <= 0.0 )
        {
            fprintf(stderr, "Invalid box dimension %d: %g\n",
                i, system_info->boxl[i] );
            exit(MD_ERR_INPUT_DATA);
        }
        system_info->boxli[i] = 1.0 / system_info->boxl[i];
    }

    // inverse mass
    for ( int i = 0; i < simulation_parameters->num_types; i++ )
    {
        if ( simulation_parameters->mass[i] <= 0.0 )
        {
            fprintf( stderr, "negative mass from input for type=%d mass=%g\n",
                i, simulation_parameters->mass[i] );
            exit(MD_ERR_INPUT_DATA);
        }
        simulation_parameters->time_over_mass[i] = 
            0.5 
            * simulation_parameters->timestep * simulation_parameters->timestep
            / simulation_parameters->mass[i];
        fprintf( stderr, "time over mass for %d is %g\n",
            i, simulation_parameters->time_over_mass[i] );
    }
    return;
}
