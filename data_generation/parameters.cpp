#include <gsl/gsl_rng.h>

#include "parameters.h"
#include "operations.h"
#include "check.h"
#include "hamiltonian.h"



void Simulation_Parameters::initialize_parameters(int number_of_occupants){
    num_occupants = number_of_occupants;
    construct_lattice(lattice);
    N = choose(num_occupants);
		b = new unsigned long long int[N]();
		table = new int[num_occupants*N]();
		ham_target = new double[2*N*N]();
		ham_initial =new double[2*N*N]();
		start_state = new double[2*N]();
    jkb_initial = new double[3]();
    jkb_target = new double[3]();
    bonds = new int[NUMBER_OF_SITES*NUMBER_OF_SITES]();

    assign_bonds(bonds, lattice);
		combinations(N, num_occupants,b,table);

    const gsl_rng_type * TT;
    gsl_rng_env_setup();
    TT = gsl_rng_default;
    rng  = gsl_rng_alloc (TT);
}



void Simulation_Parameters::initialize_target_and_initial_hamiltonian(Simulation_Parameters sim_params){
  jkb_initial[0] = g_initial;
  jkb_initial[1] = (1-g_initial)*f_initial;
  jkb_initial[2] = (1-f_initial)*(1-g_initial);

  jkb_target[0] = g_target;
  jkb_target[1] = (1-g_target)*f_target;
  jkb_target[2] = (1-f_target)*(1-g_target);


  j_initial = g_initial;
  k_initial = (1-g_initial)*f_initial;
  b_initial = (1-f_initial)*(1-g_initial);

  j_target = g_target;
  k_target = (1-g_target)*f_target;
  b_target = (1-f_target)*(1-g_target);

  construct_device_hamiltonian_uniform(sim_params, ham_initial, jkb_initial);
  construct_device_hamiltonian_uniform(sim_params, ham_target, jkb_target);

  get_ground_state(N, ham_initial, start_state);
  ground_E = get_ground_E(N, ham_target);

  if(CHECK) check_norm(start_state, N), check_commutator(N, ham_initial, ham_target);
}
