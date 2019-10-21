#ifndef __PARAMS_H_INCLUDED__
#define __PARAMS_H_INCLUDED__

#include <math.h>
#include <gsl/gsl_rng.h>


/**
    A header file containing all of the constant parameters involved in the simulation, along with the parameters class which contains all simulation
    information, and will
 */


//Physical Simulation_Parameters
const bool PERIODIC = true;
const bool UNIFORM_SITES = true;
const int DEVICE_DIMENSION = 2;
const double MAX_PARAM = 1;
const double MIN_PARAM = 0;
const double UPJ = MAX_PARAM;
const double UPK = MAX_PARAM;
const double UPB = MAX_PARAM;
const double LOWJ = MIN_PARAM;
const double LOWK = MIN_PARAM;
const double LOWB = MIN_PARAM;
const double T = 1;
const double V = 2;
const int NX = 3;
const int NY = 3;
const int NUMBER_OF_SITES = NX*NY;


//Non-Physical Simulation_Parameters
const bool MCBB = true;
const bool MCBB_DATA = true;
const bool MCBF = false;
const bool MCBF_DATA = false;
const bool ADIABATIC = false;
const bool ADIABATIC_DATA = false;
const int SEED_TOTAL= 10;
const bool DIAG = false;


//MCBF method parameters
const double DIFFERENCE_LIMIT_MC = 0.0001;
const double TAU_INIT_MC = 0.01;
const double MAX_TAU_MC = 50;
const double TAU_SCALAR_MC = 1.15;
const double MAX_CHANGE_MC_INIT = 0.5*(MAX_PARAM-MIN_PARAM);
const double ACCEPTANCE_PROB_MC = 0.8;
const double TEMP_EXP_DECAY_MC = 0.85;
const double BINARY_SEARCH_TAU_LIMIT_MC = TAU_INIT_MC/10.0;
const int RANDOM_STATES_MC = 3;
const int SWEEPS_MC = 50;
const int TOTAL_STEPS_INIT_MC =  12;
const int TEMP_DECAY_ITERATIONS_MC =  20;
const int TEMP_DECAY_LIMIT_MC = 10;
const int MAX_EVOLVE_STEPS_MC = 5*TOTAL_STEPS_INIT_MC;
const int MAX_TAU_STEPS_MC = ceil(log(MAX_TAU_MC/TAU_INIT_MC)/log(TAU_SCALAR_MC));
const int MAX_BS_STEPS_MC = ceil(log(BINARY_SEARCH_TAU_LIMIT_MC/((TAU_SCALAR_MC-1)*MAX_TAU_MC))/log(0.5));
const int ARRAY_SCALAR = 2;
const int ROW = 0;  //The change_array_mcbf variables. If 0 -> This method will not be used, if n -> use this method on every nth iteration of change_array_mcbf
const int COL = 0;
const int ALTROW = 0;
const int ALTCOL = 0;
const int ALTROW2 = 0;
const int ALTCOL2 = 0;
const int ROWK =1;
const int ROWJ =2;
const int ROWB =3;
const int COLK =0;
const int COLJ =0;
const int COLB =0;
const int SINGLE =0;
const int VARIATIONS = (ROW && 1)+(COL && 1)+(ALTROW && 1)+(ALTCOL && 1)+(ALTROW2 && 1)+(ALTCOL2 &&  1)+(SINGLE &&  1)+(ROWK && 1)+(ROWJ && 1)+(ROWB && 1)+(COLK && 1)+(COLJ && 1)+(COLB && 1);


//MCBB method parameters
const double DIFFERENCE_LIMIT_MCBB = DIFFERENCE_LIMIT_MC;
const double TAU_INIT_MCBB = TAU_INIT_MC;
const double MAX_TAU_MCBB = MAX_TAU_MC;
const double TAU_SCALAR_MCBB = TAU_SCALAR_MC;
const double CHANGE_FRACTION_MCBB = 0.1;
const double ACCEPTANCE_PROB_MCBB = ACCEPTANCE_PROB_MC;
const double TEMP_EXP_DECAY_MCBB = TEMP_EXP_DECAY_MC;
const double BINARY_SEARCH_TAU_LIMIT_MCBB = TAU_INIT_MCBB/10.0;
const int RANDOM_STATES_MCBB = RANDOM_STATES_MC;
const int NUMBER_OF_BANGS = 8;
const int SWEEPS_MCBB = 25;
const int TEMP_DECAY_ITERATIONS_MCBB = TEMP_DECAY_ITERATIONS_MC;;
const int TEMP_DECAY_LIMIT_MCBB = TEMP_DECAY_LIMIT_MC;
const int MAX_TAU_STEPS_MCBB = ceil(log(MAX_TAU_MCBB/TAU_INIT_MCBB)/log(TAU_SCALAR_MCBB));
const int MAX_BS_STEPS_MCBB = ceil(log(BINARY_SEARCH_TAU_LIMIT_MCBB/((TAU_SCALAR_MCBB-1)*MAX_TAU_MCBB))/log(0.5));


//Adiabatic method parameters
const double DIFFERENCE_LIMIT_ADIA = 0.001;
const double TAU_INIT_ADIA = 0.1;
const double MAX_TAU_ADIA = 30;
const double TAU_SCALAR_ADIA = 1.15;
const double TIME_STEP_ADIA = 1/1000.0;
const int MAX_TAU_STEPS_ADIA = ceil(MAX_TAU_ADIA/TAU_INIT_ADIA);//ceil(log(MAX_TAU_ADIA/TAU_INIT_ADIA)/log(TAU_SCALAR_ADIA));


//Debugging Help
const bool CHECK = true;
const bool PRINT = true;
const bool PRINT_COMMUTATOR = false;




class Simulation_Parameters{
public:
  double *ham_target, *ham_initial, *start_state, *state, *jkb_initial, *jkb_target, *tau_array, *best_E_array,*j_best, *k_best, *b_best;
  double ground_E, g_initial, f_initial, g_target, f_target, j_initial, k_initial, b_initial, j_target, k_target, b_target;
  double tau=0, time_step=0, temperature=0, best_E=0, difference=0, old_distance, new_distance, tau_old;
  unsigned long long int *b;
  int lattice[NX][NY], num_occupants,N,*table, *bonds, seed=0, *total_steps_cummulative_sum, destination_index, source_index, total_steps, index, best_arrays_size=0,tau_array_size=0;
  gsl_rng * rng;



  /**
      Constructs the lattice, generates the dimension of the quantum system, assigns the bonds, initializes the rng rng,
          and creates the b and table arrays (see func 'combinations' in operations.h for more info)

      @param number_of_occupants the number of occupants on the lattice
   */
  void initialize_parameters(int number_of_occupants);


  /**
      Generates the target and initial hamiltonian given the jkb initial and target values. Also grabs the ground state from
          the initial hamiltonian and the ground state from the target hamiltonian, which serve as the initial state in the
          evolution process and the target energy in the monte carlo simulations

      @param sim_params contains all of the variables for the simulation
   */
  void initialize_target_and_initial_hamiltonian(Simulation_Parameters sim_params);
};

#endif
