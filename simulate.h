/*******************************************************************************
 *
 *   Compile-time constants.  Flavor to taste
 *
 ******************************************************************************/
#define MAX_NEIGHBORS 600

#define MD_NDIMS 3

#define MD_ERR_MALLOC       1
#define MD_ERR_FOPEN        2
#define MD_ERR_INPUT_FMT    3
#define MD_ERR_INPUT_DATA   4

/******************************************************************************
 *
 *   Essential macros
 *
 ******************************************************************************/
#define COORD(i,j)      (i)*MD_NDIMS+(j)
#define LT2KT(i,j,N)      (i)*(N) - \
                            (i)*((i)+1)/2 + (j)


/******************************************************************************
 *
 *   Typedefs
 *
 ******************************************************************************/
typedef struct  atom_list               atom_list_t;
typedef struct  system_info             system_info_t;
typedef struct  sim_params              sim_params_t;
typedef struct  neighbor_list_verlet    neighbor_list_verlet_t;
typedef struct  cell_map                cell_map_t;

typedef void (*PAIR_STYLE) ( int, int, double, sim_params_t *,
                            double *, double *);
typedef void (*MULTIBODY_STYLE) ( int, int, int, double, double, 
                            sim_params_t *, double *, double *);

/******************************************************************************
 *
 *   Data Structures
 *
 ******************************************************************************/

/** @brief a collection of atoms and their identifying properties */
struct atom_list
{
    unsigned int count;     /**<@brief Number of atoms in this collection (nmol) */
    unsigned int *index;    /**<@brief Global unique identifier for each atom */
    unsigned int *type;     /**<@brief Element type for each atom */
    double      *pos,       /**<@brief Position vector on each atom */
                *vel,       /**<@brief Velocity vector on each atom */
                *accel,     /**<@brief Acceleration vector on each atom */
                *force;     /**<@brief Net force vector from interatomic potential on each atom */
    double *stress;         /**<@brief System's stress tensor, 2nd rank (3x3) */
};

/** @brief system parameters (box length, etc) */
struct system_info
{
    double boxl[MD_NDIMS];      /**<@brief Box lengths (xl, yl, zl) */
    double boxli[MD_NDIMS];     /**<@brief 1/boxl */
};

/** @brief simulation parameters such as potential parameters */
struct sim_params
{
    int num_types;          /**<@brief Number of ltypes (jtype) */
    int *ltype;             /**<@brief Maps ltypes to elements */
    double timestep;        /**<@brief integrator timestep */
    double *mass;           /**<@brief Mass of atom type (grams per atom */
    double *time_over_mass; /**<@brief timestep**2/(2*mass) */
    double **pair_param;    /**<@brief Array of pairwise potential parameters */
    PAIR_STYLE calculate_2body; /**<@brief Function to calculate 2body interactions */
    MULTIBODY_STYLE calculate_3body;    /**<@brief Function to calculate 3body interactions */
};

/** @brief the Verlet-style neighbor list */
struct neighbor_list_verlet
{
    int *nlb;
    int *nle;
    int *nlist;
};

/** @brief the link-cell method's mapping */
struct cell_map
{
    double lcell[MD_NDIMS];
    int ncell;
    int *lhead;
    int *llist;
    int *map;
};

/******************************************************************************
 *
 *   Function Prototypes
 *
 ******************************************************************************/
void init_simulation
(
    sim_params_t **simulation_parameters,
    system_info_t **system_info 
);

void simulate
(
    sim_params_t *simulation_parameters,
    system_info_t *system_info,
    atom_list_t *atoms 
);

void next_iteration
(
    sim_params_t *simulation_parameters,
    system_info_t *system_info,
    atom_list_t *atoms
);

void calculate_new_forces
(  
    sim_params_t *simulation_parameters,
    system_info_t *system_info,
    atom_list_t *atoms 
);

void simulate
(
    sim_params_t *simulation_parameters,
    system_info_t *system_info,
    atom_list_t *atoms 
);

void pair_lj
( 
    int itype, 
    int jtype, 
    double rij2, 
    sim_params_t *simulation_parameters,
    double *potential_mag, 
    double *force_mag 
);

neighbor_list_verlet_t *construct_neighbor_list
(  
    system_info_t *system_info,
    atom_list_t *atoms 
);

void malloc_error( void );

atom_list_t *get_starting_config
(
    sim_params_t *simulation_parameters,
    system_info_t *system_info 
);

atom_list_t *get_starting_config_xyz
(
    sim_params_t *simulation_parameters,
    system_info_t *system_info 
);

cell_map_t *create_link_cell
( 
    system_info_t *system_info, 
    atom_list_t *atoms 
);

void init_constants
(
    sim_params_t *simulation_parameters,
    system_info_t *system_info 
);

void modify_starting_config( system_info_t *system_info, atom_list_t *atoms );
void get_input_params( sim_params_t *simulation_parameters );
void dump_table( sim_params_t *simulation_parameters );
