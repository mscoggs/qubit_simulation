#include <gsl/gsl_rng.h>
#include <cstring>
#include <algorithm>    // std::min

#include "check.h"
#include "export_data.h"
#include "hamiltonian.h"
#include "linear_algebra.h"
#include "mcdb.h"
#include "operations.h"
#include "parameters.h"
#include "print.h"




void mcdb_method(Simulation_Parameters& sim_params){
	int i;
	sim_params.init_mcdb_params();

	if(check_commutator(sim_params.N, sim_params.ham_initial, sim_params.ham_target) || sim_params.initial_E -sim_params.ground_E < 0.001){
		sim_params.tau = 0.0, sim_params.new_distance = 0.0, sim_params.best_E = 0.0;
		if(PRINT) print_mcdb_info(sim_params);
		if(MCDB_DATA) save_mcdb_data_fixed_tau(sim_params);
		sim_params.clear_mcdb_params();
		return;
	}

	pre_exponentiate(sim_params);

	while(sim_params.tau<MAX_TAU_MCDB){

		gsl_rng_set(sim_params.rng, 1);
		calc_initial_temp_mcdb(sim_params);

		for(sim_params.seed=1; sim_params.seed<NUM_SEEDS+1; sim_params.seed++){
			gsl_rng_set(sim_params.rng, sim_params.seed);
			mcdb_simulation(sim_params);

			sim_params.E_array_fixed_tau[sim_params.seed-1] = sim_params.best_E;
			copy_arrays_mcdb(sim_params,sim_params.j_best_fixed_tau,sim_params.k_best_fixed_tau,sim_params.b_best_fixed_tau,sim_params.j_best,  sim_params.k_best,  sim_params.b_best,  (sim_params.seed-1)*sim_params.total_steps, 0);
		}
		for(sim_params.seed=1; sim_params.seed<NUM_SEEDS+1; sim_params.seed++) if(sim_params.E_array_fixed_tau[sim_params.seed-1] < sim_params.best_E) sim_params.best_E = sim_params.E_array_fixed_tau[sim_params.seed-1];

		if(!update_distances(sim_params)) continue;

		if(PRINT) print_mc_results(sim_params);
    if(MCDB_DATA) save_mcdb_data_fixed_tau(sim_params);

		if(sim_params.new_distance < DISTANCE_LIMIT_MCDB) break;
		else{
      calc_tau_mcdb(sim_params);
      sim_params.total_steps = (int)floor(sim_params.tau/sim_params.time_step);
		}
	}
	sim_params.clear_mcdb_params();

}




void mcdb_simulation(Simulation_Parameters& sim_params){
	int i,j, random_time_index, proposal_accepted,proposal_count, poor_acceptance_streak;
	double *j_array, *k_array,*b_array,*j_temp, *k_temp, *b_temp, acceptance_rate, old_E, new_E, change;


	j_array          = new double[sim_params.total_steps]();
	k_array          = new double[sim_params.total_steps]();
	b_array          = new double[sim_params.total_steps]();
	j_temp           = new double[sim_params.total_steps]();
	k_temp           = new double[sim_params.total_steps]();
	b_temp           = new double[sim_params.total_steps]();
	sim_params.state = new double[2*sim_params.N]();


	std::memcpy(sim_params.state,sim_params.start_state, 2*sim_params.N*sizeof(double));

	init_arrays_mcdb(sim_params, sim_params.j_best,sim_params.k_best,sim_params.b_best);
	copy_arrays_mcdb(sim_params, j_array, k_array, b_array,sim_params.j_best,sim_params.k_best,sim_params.b_best,0,0);//a temporary array, used in the undoing of the changes

	evolve_mcdb(sim_params, j_array,k_array,b_array);
	sim_params.best_E = cost(sim_params.N, sim_params.state, sim_params.ham_target);
	old_E = sim_params.best_E;

	if(PRINT) print_mcdb_info(sim_params);

	for (i=0;i<TEMP_DECAY_ITERATIONS_MCDB;i++){

		proposal_accepted = 0, proposal_count = 0;

		for (j=0; j<SWEEPS_MCDB*sim_params.total_steps*sim_params.sweeps_multiplier;j++){

			copy_arrays_mcdb(sim_params, j_temp, k_temp,b_temp,j_array, k_array, b_array,0,0);//a temporary array, used in the undoing of the changes
			random_time_index = (int)floor(get_random_double(0, sim_params.total_steps, sim_params.rng));
			change_array_mcdb(j_array, k_array, b_array,random_time_index, j);

			std::memcpy(sim_params.state,sim_params.start_state, 2*sim_params.N*sizeof(double));//resetting state
			evolve_mcdb(sim_params, j_array, k_array, b_array);
			new_E = cost(sim_params.N,sim_params.state, sim_params.ham_target);

			if (new_E<sim_params.best_E) sim_params.best_E=new_E, old_E=new_E,  copy_arrays_mcdb(sim_params, sim_params.j_best, sim_params.k_best, sim_params.b_best,j_array, k_array, b_array, 0, 0), poor_acceptance_streak = 0;
			//else if (new_E<=old_E) old_E=new_E;
			else if (get_random_double(0,1,sim_params.rng)<exp(-(new_E-old_E)/(sim_params.temperature))) old_E=new_E, proposal_accepted++, proposal_count++;
			else copy_arrays_mcdb(sim_params, j_array, k_array, b_array,  j_temp, k_temp, b_temp, 0,0), proposal_count++;//undoing the change
		}

		acceptance_rate = (double)proposal_accepted/proposal_count;
		if(acceptance_rate<0.1) poor_acceptance_streak++;
		else poor_acceptance_streak = 0;

		if(PRINT) printf("           Best Expectation:   %3.6f  ||  Acceptance Rate: %3.4f (%i/%i)\n",sim_params.best_E,acceptance_rate,proposal_accepted, proposal_count);

		if(poor_acceptance_streak>TEMP_DECAY_LIMIT_MCDB){
			printf("NO MC PROGRESS FOR %i TEMP DECAY LIMIT MCDB, TERMINATING\n", TEMP_DECAY_LIMIT_MCDB);
			break;
		}

		sim_params.temperature=sim_params.temperature*TEMP_EXP_DECAY_MCDB;
	}
	delete[] k_array, delete[] j_array, delete[] b_array, delete[] k_temp, delete[] j_temp, delete[] b_temp, delete[] sim_params.state;
}




void evolve_mcdb(Simulation_Parameters& sim_params, double *j_array, double *k_array, double *b_array){
	int i;

	for(i=0; i<sim_params.total_steps; i++){
		//print_state(sim_params.state, sim_params.N);

    if(j_array[i] == 1.0 && k_array[i] == 1.0)	matrix_vector_mult(sim_params.e11,sim_params.state, sim_params.N);
    else if(j_array[i] == 0.0 && k_array[i] == 1.0)	matrix_vector_mult(sim_params.e01,sim_params.state, sim_params.N);
    else if(j_array[i] == 1.0 && k_array[i] == 0.0)	matrix_vector_mult(sim_params.e10,sim_params.state, sim_params.N);
    else if(j_array[i] == 0.0 && k_array[i] == 0.0) continue;
    else printf("NOT EVOLVING WHAT'S GOING ON?!?!?\n\n\n\n\n"), exit(0);

  if(CHECK) check_norm(sim_params.state, sim_params.N);
  }
}




void calc_initial_temp_mcdb(Simulation_Parameters& sim_params){
	int i, j, random_time_index, start_state_index, count=0;
	double *state_random, *j_temp, *k_temp, *b_temp, sum=0, change, old_E, new_E;

	j_temp           = new double[sim_params.total_steps]();
	k_temp           = new double[sim_params.total_steps]();
	b_temp           = new double[sim_params.total_steps]();
	sim_params.state = new double[2*sim_params.N]();
	state_random     = new double[2*sim_params.N]();

	if(PRINT) printf("\n\n\n\n...Calculating initial temperature based on %i random starting states...\n", RANDOM_STATES_MCDB);

	for (j=0;j<RANDOM_STATES_MCDB;j++){
		for (i=0; i<sim_params.N*2;i++) state_random[i] =0.0;
		start_state_index = floor(get_random_double(0,sim_params.N,sim_params.rng));
		state_random[start_state_index*2] = 1;
		std::memcpy(sim_params.state,state_random, 2*sim_params.N*sizeof(double));

		init_arrays_mcdb(sim_params, j_temp, k_temp, b_temp);
		evolve_mcdb(sim_params,j_temp, k_temp, b_temp);
		old_E = cost(sim_params.N, sim_params.state, sim_params.ham_target);

		for (i=0; i<SWEEPS_MCDB*sim_params.total_steps;i++){
			random_time_index = (int)floor(get_random_double(0, sim_params.total_steps, sim_params.rng));
      change_array_mcdb(j_temp, k_temp, b_temp,random_time_index, i);

			std::memcpy(sim_params.state,state_random, 2*sim_params.N*sizeof(double));//resetting state
			evolve_mcdb(sim_params,j_temp, k_temp, b_temp);
			new_E = cost(sim_params.N, sim_params.state, sim_params.ham_target);

			if (new_E>=old_E) sum += (new_E-old_E), count++;
			old_E=new_E;
		}
	}
	sim_params.temperature = -(sum/(count*log(ACCEPTANCE_PROB_MC)));

	delete[] j_temp, delete[] k_temp, delete[] b_temp, delete[] sim_params.state, delete[] state_random;
}



void calc_tau_mcdb(Simulation_Parameters& sim_params){
	double tau_scalar;


	if(sim_params.new_distance < 0.2) tau_scalar = TAU_SCALAR_MCDB_TINY;
	else if((abs(sim_params.old_distance - sim_params.new_distance) < .1) and (sim_params.new_distance > 0.2 )) tau_scalar = TAU_SCALAR_MCDB_BIG;
	else tau_scalar = TAU_SCALAR_MCDB;

	sim_params.tau = sim_params.tau*tau_scalar;
}




void init_arrays_mcdb(Simulation_Parameters& sim_params, double *j_array,double *k_array,double *b_array){
	int i;
	for (i=0; i<sim_params.total_steps;i++){
		      j_array[i] = fmod(i+2.0, 2);//round(get_random_double(0,1.0,sim_params.rng));
	       	k_array[i] = fmod(i+1.0, 2);//round(get_random_double(0,1.0,sim_params.rng));
	       	b_array[i] = 1.0;//round(get_random_double(0,1.0,sim_params.rng));
	}
}




void copy_arrays_mcdb(Simulation_Parameters& sim_params, double* j_to,  double* k_to, double* b_to, double *j_from, double *k_from, double* b_from, int destination_index, int source_index){


	std::memcpy(&j_to[destination_index], &j_from[source_index], sizeof(double)*sim_params.total_steps);
	std::memcpy(&k_to[destination_index], &k_from[source_index], sizeof(double)*sim_params.total_steps);
	std::memcpy(&b_to[destination_index], &b_from[source_index], sizeof(double)*sim_params.total_steps);
}




void change_array_mcdb(double *j_array, double *k_array, double *b_array, int random_time_index, int i){
	double *pointer, *pointer2;

	if(i%2 == 0) pointer = j_array, pointer2 = k_array;
	if(i%2 == 1) pointer = k_array, pointer2 = j_array;
	//if(i%3 == 2) pointer = b_array;

  *(pointer+random_time_index) =  fmod(*(pointer+random_time_index) + 1.0, 2.0);
	if(*(pointer+random_time_index)+*(pointer2+random_time_index)> 0.01) return;
	else *(pointer+random_time_index) += 1;
}



void pre_exponentiate(Simulation_Parameters& sim_params){
	int j;
	double *hamiltonian,*ham_t_i,*jkb_vals;

	hamiltonian  = new double[2*sim_params.N*sim_params.N]();
	ham_t_i      = new double[2*sim_params.N*sim_params.N]();
	jkb_vals     = new double[3]();
  for (j=0; j<sim_params.N*sim_params.N*2; j++) sim_params.e11[j]=0, sim_params.e01[j]=0,sim_params.e10[j]=0;

	jkb_vals[0] = 1, jkb_vals[1] = 1, jkb_vals[2] = 0;
	construct_device_hamiltonian_uniform(sim_params, hamiltonian, jkb_vals);
	for (j=0; j<sim_params.N*sim_params.N; j++) ham_t_i[2*j+1] = (hamiltonian[2*j]*-sim_params.time_step), ham_t_i[2*j] = (hamiltonian[2*j+1]*sim_params.time_step);//multiplying by -i*dt for the Pade approximation
	exp_complex_double_matrix_pade(sim_params.N, ham_t_i, sim_params.e11);

	jkb_vals[0] = 0, jkb_vals[1] = 1, jkb_vals[2] = 0;
	construct_device_hamiltonian_uniform(sim_params, hamiltonian, jkb_vals);
	for (j=0; j<sim_params.N*sim_params.N; j++) ham_t_i[2*j+1] = (hamiltonian[2*j]*-sim_params.time_step), ham_t_i[2*j] = (hamiltonian[2*j+1]*sim_params.time_step);//multiplying by -i*dt for the Pade approximation
	exp_complex_double_matrix_pade(sim_params.N, ham_t_i, sim_params.e01);

	jkb_vals[0] = 1.0, jkb_vals[1] = 0, jkb_vals[2] = 0;
	construct_device_hamiltonian_uniform(sim_params, hamiltonian, jkb_vals);
	for (j=0; j<sim_params.N*sim_params.N; j++) ham_t_i[2*j+1] = (hamiltonian[2*j]*-sim_params.time_step), ham_t_i[2*j] = (hamiltonian[2*j+1]*sim_params.time_step);//multiplying by -i*dt for the Pade approximation
	exp_complex_double_matrix_pade(sim_params.N, ham_t_i, sim_params.e10);
	if(CHECK){
		check_unitary(sim_params.e11, sim_params.N);
		check_unitary(sim_params.e10, sim_params.N);
		check_unitary(sim_params.e01, sim_params.N);
	}

	delete[] jkb_vals, delete[] ham_t_i, delete[] hamiltonian;
}
