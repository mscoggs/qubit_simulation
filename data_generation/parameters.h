#ifndef __PARAMS_H_INCLUDED__
#define __PARAMS_H_INCLUDED__

#include <ctime>
#include <math.h>
#include <gsl/gsl_rng.h>




/**
    A header file containing all of the constant parameters involved in the simulation, along with the parameters class which contains all simulation variables
 */



/*SWITCHES FOR EACH TYPE OF SIMULATION METHOD*/
 const bool MCDB      = true;
 const bool MCBB      = false;
 const bool MCBF      = false;
 const bool ADIA      = false;
 const bool SAVE_DATA = true;


/*LATTICE PARAMETERS*/
const bool   USE_ENERGY_DISTANCE   = false;
const bool   PERIODIC              = false;
const bool   UNIFORM_SITES         = true;
const double MAX_PARAM             = 1.0;
const double MIN_PARAM             = 0.0;
const int    DEVICE_DIMENSION      = 2;
const int    NX                    = 3;
const int    NY                    = NX;
const int    NUMBER_OF_SITES       = NX*NY;

/*SIMULATION PARAMETERS*/
const bool   DIAG                  = true;
const int    NUM_SEEDS             = 1;
const double DISTANCE_LIMIT        = 0.01;
const double TAU_INIT              = 0.1;
const double MAX_TAU               = 3;
const double TAU_SCALAR            = 1.15;
const double TAU_SCALAR_TINY       = 1.05;
const double TAU_SCALAR_BIG        = 1.3;
const double ACCEPTANCE_PROB       = 0.65;
const double TEMP_EXP_DECAY        = 0.90;
const double MIN_TEMP_FRACTION     = 0.01;
const int    TEMP_DECAY_ITERATIONS = ceil(log(MIN_TEMP_FRACTION)/log(TEMP_EXP_DECAY)); //30 for our given values
const int    ZERO_TEMP_ITERATIONS  = 15;
const int    RANDOM_STATES         = 4;
const double INIT_OVERLAP_LIMIT    = 0.99999;


/*MCBB METHOD PARAMETERS*/
const double MAX_CHANGE_FRACTION_MCBB = 0.9;
const double MIN_CHANGE_FRACTION_MCBB = 0.02;
const int    NUMBER_OF_BANGS_MCBB     = 6;
const int    SWEEPS_MCBB              = 90;

/*MCDB METHOD PARAMETERS*/
const int    MAX_STEPS_MCDB    = 128;
const int    MIN_STEPS_MCDB    = 4; //MAKE SURE THIS IS LESS THAN OR EQUAL TO THE NUMBER OF BANGS
const int    NUMBER_OF_BANGS_MCDB     = MIN_STEPS_MCDB;
const int    TOTAL_STEP_CHANGES= (int)round((log2(MAX_STEPS_MCDB))) + 1;
const int    SWEEPS_MCDB       = 50;
const double STEPS_CRUNCH_MCDB = 1.0;

/*MCBF METHOD PARAMETERS*/
const double MAX_CHANGE_MCBF = 0.6*(MAX_PARAM-MIN_PARAM);
const double MIN_CHANGE_MCBF = 0.05*(MAX_PARAM-MIN_PARAM);
const int    SWEEPS_MCBF            = 50;
const int    TOTAL_STEPS_INIT_MCBF  = 5;
const int    MAX_EVOLVE_STEPS_MCBF  = 4*TOTAL_STEPS_INIT_MCBF;
const int    ARRAY_SCALAR         = 2;
const int ROW        = 0;  //The change_array_mcbf variables. If 0 -> This method will not be used, if n -> use this method on every nth iteration of change_array_mcbf
const int COL        = 0;
const int ALTROW     = 0;
const int ALTCOL     = 0;
const int ALTROW2    = 0;
const int ALTCOL2    = 0;
const int ROWK       = 1; //Rows correspond to all sites for a fixed time, columns to all times for a fixed site.
const int ROWJ       = 2;
const int ROWB       = 0;
const int COLK       = 0;
const int COLJ       = 0;
const int COLB       = 0;
const int SINGLE     = 0;
const int VARIATIONS = (ROW && 1)+(COL && 1)+(ALTROW && 1)+(ALTCOL && 1)+(ALTROW2 && 1)+(ALTCOL2 &&  1)+(SINGLE &&  1)+(ROWK && 1)+(ROWJ && 1)+(ROWB && 1)+(COLK && 1)+(COLJ && 1)+(COLB && 1);

/*ADIABATIC METHOD PARAMETERS*/
const double MAX_TAU_ADIA    = 200;
const double TAU_SCALAR_ADIA = 1.5;
const int TOTAL_STEPS_ADIA   = 50;



//Tests and debugging help
const bool CHECK            = true;
const bool PRINT            = true;
const bool PRINT_COMMUTATOR = false;




class Simulation_Parameters{
public:
	std::clock_t start;
	double *ham_target, *ham_initial, *init_state, *target_state, *state, *jkb_initial, *jkb_target, *best_mc_result_fixed_tau, *evolved_state_fixed_tau, *best_evolved_state, *j_best_fixed_tau, *k_best_fixed_tau, *b_best_fixed_tau, *j_best, *k_best, *b_best, *e11, *e01, *e10, *j_best_scaled, *k_best_scaled, *b_best_scaled, *saved_states, *P_11, *P_11_inv, *D_11, *P_10, *P_10_inv, *D_10, *P_01, *P_01_inv, *D_01;
	double ground_E, j_initial, k_initial, b_initial, j_target, k_target, b_target, tau, time_step, temperature,initial_temperature, best_mc_result,initial_E, old_distance, new_distance, temp_distance, init_target_dot_squared, evolved_target_dot_squared, duration;
	unsigned long long int *b;
	int lattice[NX][NY], num_occupants,N,*table, *bonds, seed, total_steps, max_steps_mcdb, sweeps_multiplier, total_sweeps, j_bangs, k_bangs, max_jumps;
	gsl_rng * rng;



	 /**
	 Constructs the lattice, generates the dimension of the quantum system, assigns the bonds, initializes the rng rng,
	 and creates the b and table arrays (see func 'combinations' in operations.h for more info)

	 @param number_of_occupants occupants on the lattice
	 @param sim_params contains all of the variables for the simulation
	 */
	void initialize_lattice(int number_of_occupants, Simulation_Parameters& sim_params);



	/**
	Generates the target and initial hamiltonian given the jkb initial and target values. Also grabs the ground state from
	the initial hamiltonian and the ground state from the target hamiltonian, which serve as the initial state in the
	evolution process and the target energy in the monte carlo simulations

	@param ji j_initial
	@param ki k_initial
	@param jt j_target
	@param kt k_target
	@param sim_params contains all of the variables for the simulation
	*/
	void initialize_hamiltonians(double ji, double ki, double jt, double kt, Simulation_Parameters& sim_params);


	/**
		Initializing/clearing the arrays for the simulations
	*/
	void init_mcbb_params(Simulation_Parameters& sim_params);
	void clear_mcbb_params();
	void init_mcdb_params();
	void clear_mcdb_params();
	void init_mcbf_params();
	void clear_mcbf_params();
};




#endif
