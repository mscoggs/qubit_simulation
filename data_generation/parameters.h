#ifndef __PARAMS_H_INCLUDED__
#define __PARAMS_H_INCLUDED__

#include <math.h>
#include <gsl/gsl_rng.h>




/**
    A header file containing all of the constant parameters involved in the simulation, along with the parameters class which contains all simulation variables
 */




//Physical Simulation_Parameters
const bool PERIODIC = true;
const bool UNIFORM_SITES = true;
const int DEVICE_DIMENSION = 2;
const double MAX_PARAM = 1.0;
const double MIN_PARAM = 0.0;
const double UPJ = MAX_PARAM;
const double UPK = MAX_PARAM;
const double UPB = MAX_PARAM;
const double LOWJ = MIN_PARAM;
const double LOWK = MIN_PARAM;
const double LOWB = MIN_PARAM;
const double T = 1;
const double V = 2;
const int NX = 2;
const int NY = NX;
const int NUMBER_OF_SITES = NX*NY;



//Non-Physical Simulation_Parameters
const bool MCBF = false;
const bool MCBF_DATA = false;
const bool MCBB = true;
const bool MCBB_DATA = true;
const bool MCDB = false;
const bool MCDB_DATA = false;
const bool ADIA = false;
const bool ADIA_DATA= false;
const int  NUM_SEEDS = 4;
const bool DIAG = false;



//MCBF method parameters
const double DISTANCE_LIMIT_MCBF = 0.02;
const double TAU_INIT_MCBF = 0.01;
const double MAX_TAU_MCBF = 0.7;
const double TAU_SCALAR_MCBF = 1.3;
const double TAU_SCALAR_MCBF_TINY = 1.1;
const double TAU_SCALAR_MCBF_BIG = 1.8;
const double MAX_CHANGE_MCBF_INIT = 0.5*(MAX_PARAM-MIN_PARAM);
const double MIN_CHANGE_MCBF_INIT = 0.1*(MAX_PARAM-MIN_PARAM);
const double ACCEPTANCE_PROB_MC = 0.80;
const double TEMP_EXP_DECAY_MC = 0.80;
const double BINARY_SEARCH_TAU_LIMIT_MCBF = TAU_INIT_MCBF/100.0;
const int RANDOM_STATES_MC = 3;
const int SWEEPS_MC = 50;
const int TOTAL_STEPS_INIT_MC =  5;
const int TEMP_DECAY_ITERATIONS_MC = 20;
const int TEMP_DECAY_LIMIT_MC = 7;
const int MAX_EVOLVE_STEPS_MC = 4*TOTAL_STEPS_INIT_MC;
const int MAX_TAU_STEPS_MCBF = ceil(log(MAX_TAU_MCBF/TAU_INIT_MCBF)/log(TAU_SCALAR_MCBF));
const int MAX_BS_STEPS_MCBF = ceil(log(BINARY_SEARCH_TAU_LIMIT_MCBF/((TAU_SCALAR_MCBF-1)*MAX_TAU_MCBF))/log(0.5));
const int ARRAY_SCALAR = 2;
const int ROW = 0;  //The change_array_mcbf variables. If 0 -> This method will not be used, if n -> use this method on every nth iteration of change_array_mcbf
const int COL = 0;
const int ALTROW = 0;
const int ALTCOL = 0;
const int ALTROW2 = 0;
const int ALTCOL2 = 0;
const int ROWK =1;
const int ROWJ =2;
const int ROWB =0;
const int COLK =0;
const int COLJ =0;
const int COLB =0;
const int SINGLE =0;
const int VARIATIONS = (ROW && 1)+(COL && 1)+(ALTROW && 1)+(ALTCOL && 1)+(ALTROW2 && 1)+(ALTCOL2 &&  1)+(SINGLE &&  1)+(ROWK && 1)+(ROWJ && 1)+(ROWB && 1)+(COLK && 1)+(COLJ && 1)+(COLB && 1);



//MCBB method parameters
const double DISTANCE_LIMIT_MCBB = DISTANCE_LIMIT_MCBF;
const double TAU_INIT_MCBB = TAU_INIT_MCBF;
const double MAX_TAU_MCBB = MAX_TAU_MCBF;
const double TAU_SCALAR_MCBB = TAU_SCALAR_MCBF;
const double TAU_SCALAR_MCBB_TINY = TAU_SCALAR_MCBF_TINY;
const double TAU_SCALAR_MCBB_BIG = TAU_SCALAR_MCBF_BIG;
const double MAX_CHANGE_FRACTION_MCBB = 0.8;
const double MIN_CHANGE_FRACTION_MCBB = 0.01;
const double ACCEPTANCE_PROB_MCBB = ACCEPTANCE_PROB_MC;
const double TEMP_EXP_DECAY_MCBB = TEMP_EXP_DECAY_MC;
const double BINARY_SEARCH_TAU_LIMIT_MCBB =BINARY_SEARCH_TAU_LIMIT_MCBF;
const int RANDOM_STATES_MCBB = RANDOM_STATES_MC;
const int NUMBER_OF_BANGS = 3;
const int SWEEPS_MCBB = 50;
const int TEMP_DECAY_ITERATIONS_MCBB = TEMP_DECAY_ITERATIONS_MC;;
const int TEMP_DECAY_LIMIT_MCBB = TEMP_DECAY_LIMIT_MC;
const int MAX_TAU_STEPS_MCBB = ceil(log(MAX_TAU_MCBB/TAU_INIT_MCBB)/log(TAU_SCALAR_MCBB));
const int MAX_BS_STEPS_MCBB = ceil(log(BINARY_SEARCH_TAU_LIMIT_MCBB/((TAU_SCALAR_MCBB-1)*MAX_TAU_MCBB))/log(0.5));



//MCDB method parameters
const double DISTANCE_LIMIT_MCDB = DISTANCE_LIMIT_MCBF;
const double TAU_INIT_MCDB = TAU_INIT_MCBF;
const double MAX_TAU_MCDB = MAX_TAU_MCBF;
const double TAU_SCALAR_MCDB = TAU_SCALAR_MCBF;
const double TAU_SCALAR_MCDB_TINY = TAU_SCALAR_MCBF_TINY;
const double TAU_SCALAR_MCDB_BIG = TAU_SCALAR_MCBF_BIG;
const double ACCEPTANCE_PROB_MCDB = ACCEPTANCE_PROB_MC;
const double TEMP_EXP_DECAY_MCDB = TEMP_EXP_DECAY_MC;
const double BINARY_SEARCH_TAU_LIMIT_MCDB =BINARY_SEARCH_TAU_LIMIT_MCBF;
const int RANDOM_STATES_MCDB = RANDOM_STATES_MC;
const double TIME_STEP_MCDB = 0.005;
const int SWEEPS_MCDB = 100;
const int TEMP_DECAY_ITERATIONS_MCDB = TEMP_DECAY_ITERATIONS_MC;;
const int TEMP_DECAY_LIMIT_MCDB = TEMP_DECAY_LIMIT_MC;
const int MAX_TAU_STEPS_MCDB = ceil(log(MAX_TAU_MCDB/TAU_INIT_MCDB)/log(TAU_SCALAR_MCDB));
const int MAX_BS_STEPS_MCDB = ceil(log(BINARY_SEARCH_TAU_LIMIT_MCDB/((TAU_SCALAR_MCDB-1)*MAX_TAU_MCDB))/log(0.5));





//Adiabatic method parameters
const double DISTANCE_LIMIT_ADIA = DISTANCE_LIMIT_MCBF;
const double TAU_INIT_ADIA = TAU_INIT_MCBF;
const double MAX_TAU_ADIA = MAX_TAU_MCBF;
const double TAU_SCALAR_ADIA = 1.15;
const double TIME_STEP_ADIA = 1/1000.0;
const int MAX_TAU_STEPS_ADIA = ceil(MAX_TAU_ADIA/TAU_INIT_ADIA);//ceil(log(MAX_TAU_ADIA/TAU_INIT_ADIA)/log(TAU_SCALAR_ADIA));



//Tests and debugging help
const bool CHECK = true;
const bool PRINT = true;
const bool PRINT_COMMUTATOR = false;




class Simulation_Parameters{
public:
	double *ham_target, *ham_initial, *start_state, *state, *jkb_initial, *jkb_target, *E_array_fixed_tau, *j_best_fixed_tau, *k_best_fixed_tau, *b_best_fixed_tau, *j_best, *k_best, *b_best, *e11, *e01, *e10;
	double ground_E, j_initial, k_initial, b_initial, j_target, k_target, b_target, tau, time_step, temperature, best_E,initial_E, old_distance, new_distance, temp_distance;
	unsigned long long int *b;
	int lattice[NX][NY], num_occupants,N,*table, *bonds, seed, total_steps, temp_iteration, max_steps_mcdb, sweeps_multiplier;
	gsl_rng * rng;




	 /**
	 Constructs the lattice, generates the dimension of the quantum system, assigns the bonds, initializes the rng rng,
	 and creates the b and table arrays (see func 'combinations' in operations.h for more info)

		Generates the target and initial hamiltonian given the jkb initial and target values. Also grabs the ground state from
		the initial hamiltonian and the ground state from the target hamiltonian, which serve as the initial state in the
		evolution process and the target energy in the monte carlo simulations

		@param sim_params contains all of the variables for the simulation
	 */
	void initialize_simluation(int number_of_occupants, double ji, double ki, double jt, double kt, Simulation_Parameters& sim_params);




	/**
		Initializing/clearing the arrays for the simulations
	*/
	void init_mcbb_params();
	void clear_mcbb_params();
	void init_mcdb_params();
	void clear_mcdb_params();
	void init_mcbf_params();
	void clear_mcbf_params();
};




#endif
