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
		if(SAVE_DATA) save_mcdb_data_fixed_tau(sim_params);
		sim_params.clear_mcdb_params();
		return;
	}


	while(sim_params.tau<MAX_TAU){

		sim_params.total_steps = MIN_STEPS_MCDB;
		sim_params.time_step = sim_params.tau/sim_params.total_steps;

		while(sim_params.total_steps <= MAX_STEPS_MCDB){

			pre_exponentiate(sim_params);
			gsl_rng_set(sim_params.rng, 0);
			get_sweeps_mcdb(sim_params);
			calc_initial_temp_mcdb(sim_params);

			for(sim_params.seed=1; sim_params.seed<NUM_SEEDS+1; sim_params.seed++){
				gsl_rng_set(sim_params.rng, sim_params.seed);
				mcdb_simulation(sim_params);

				sim_params.E_array_fixed_tau[sim_params.seed-1] = sim_params.best_E;
				copy_arrays_mcdb(sim_params,sim_params.j_best_fixed_tau,sim_params.k_best_fixed_tau,sim_params.b_best_fixed_tau,sim_params.j_best,  sim_params.k_best,  sim_params.b_best,  (sim_params.seed-1)*sim_params.total_steps, 0);
			}

			for(sim_params.seed=1; sim_params.seed<NUM_SEEDS+1; sim_params.seed++){
				if(sim_params.E_array_fixed_tau[sim_params.seed-1] <= sim_params.best_E){
				       	sim_params.best_E = sim_params.E_array_fixed_tau[sim_params.seed-1];
					copy_arrays_mcdb(sim_params,sim_params.j_best_scaled,  sim_params.k_best_scaled,  sim_params.b_best_scaled,sim_params.j_best_fixed_tau,sim_params.k_best_fixed_tau,sim_params.b_best_fixed_tau,  0,(sim_params.seed-1)*sim_params.total_steps);
				}
			}

			if(sim_params.total_steps == MAX_STEPS_MCDB) break;

			if(sim_params.new_distance > 0.2 and sim_params.total_steps >= MIN_STEPS_MCDB*4 and sim_params.total_steps <= MAX_STEPS_MCDB/2) break;
			sim_params.total_steps = sim_params.total_steps*2;
			sim_params.time_step = sim_params.tau/sim_params.total_steps;
			if(sim_params.total_steps <= MAX_STEPS_MCDB) scale_best_arrays_mcdb(sim_params, sim_params.j_best_scaled,sim_params.k_best_scaled,sim_params.b_best_scaled);

		}

		if(!update_distances(sim_params)) continue;
		if(PRINT) print_mc_results(sim_params);
		if(SAVE_DATA) {
			sim_params.duration = (std::clock() - sim_params.start)/(double) CLOCKS_PER_SEC;
			save_mcdb_data_fixed_tau(sim_params);
		}
		if(sim_params.old_distance < sim_params.new_distance) break;
		if(sim_params.new_distance < DISTANCE_LIMIT) break;
		calc_tau(sim_params);
	}
	sim_params.clear_mcdb_params();
}




void mcdb_simulation(Simulation_Parameters& sim_params){
	int i,j, k,random_time_index, proposal_accepted,proposal_count, poor_acceptance_streak;
	double *j_array, *k_array,*b_array,*j_temp, *k_temp, *b_temp,*temp_saved_states, acceptance_rate, old_E, new_E, change;


	j_array          = new double[sim_params.total_steps]();
	k_array          = new double[sim_params.total_steps]();
	b_array          = new double[sim_params.total_steps]();
	j_temp           = new double[sim_params.total_steps]();
	k_temp           = new double[sim_params.total_steps]();
	b_temp           = new double[sim_params.total_steps]();
	sim_params.state = new double[2*sim_params.N]();
	temp_saved_states = new double[2*sim_params.N*sim_params.total_steps]();



	if(sim_params.total_steps == MIN_STEPS_MCDB) init_arrays_mcdb(sim_params, sim_params.j_best,sim_params.k_best,sim_params.b_best);
	else copy_arrays_mcdb(sim_params,sim_params.j_best,sim_params.k_best,sim_params.b_best, sim_params.j_best_scaled,sim_params.k_best_scaled,sim_params.b_best_scaled,0,0);
	copy_arrays_mcdb(sim_params, j_array, k_array, b_array,sim_params.j_best,sim_params.k_best,sim_params.b_best,0,0);//a temporary array, used in the undoing of the changes

	std::memcpy(sim_params.state,sim_params.start_state, 2*sim_params.N*sizeof(double));
	evolve_mcdb(sim_params, j_array,k_array,b_array,0);
	sim_params.best_E = cost(sim_params.N, sim_params.state, sim_params.ham_target);
	old_E = sim_params.best_E;

	if(PRINT) print_mcdb_info(sim_params);


	sim_params.temperature = sim_params.initial_temperature;

	for (i=0;i<TEMP_DECAY_ITERATIONS;i++){

		proposal_accepted = 0, proposal_count = 0;

		for (j=0; j<sim_params.total_sweeps*sim_params.sweeps_multiplier;j++){

			copy_arrays_mcdb(sim_params, j_temp, k_temp,b_temp,j_array, k_array, b_array,0,0);//a temporary array, used in the undoing of the changes
			std::memcpy(temp_saved_states,sim_params.saved_states, sim_params.total_steps*2*sim_params.N*sizeof(double));

			random_time_index = change_array_mcdb(sim_params, j_array, k_array, b_array, j);

			if(j==0) random_time_index = 0;
			for(k=0; k<2*sim_params.N; k++) sim_params.state[k] = sim_params.saved_states[(random_time_index)*2*sim_params.N +k];
			evolve_mcdb(sim_params, j_array, k_array, b_array, random_time_index);

			new_E = cost(sim_params.N,sim_params.state, sim_params.ham_target);

			if (new_E<sim_params.best_E) sim_params.best_E=new_E, old_E=new_E,  copy_arrays_mcdb(sim_params, sim_params.j_best, sim_params.k_best, sim_params.b_best,j_array, k_array, b_array, 0, 0), poor_acceptance_streak = 0;
			else if (new_E<=old_E) old_E=new_E;
			else if (get_random_double(0,1,sim_params.rng)<exp(-(new_E-old_E)/(sim_params.temperature))) old_E=new_E, proposal_accepted++, proposal_count++;
			else copy_arrays_mcdb(sim_params, j_array, k_array, b_array,  j_temp, k_temp, b_temp, 0,0), std::memcpy(sim_params.saved_states, temp_saved_states,sim_params.total_steps*2*sim_params.N*sizeof(double)), proposal_count++;//undoing the change
		}

		acceptance_rate = (double)proposal_accepted/proposal_count;
		if(acceptance_rate<0.1) poor_acceptance_streak++;
		else poor_acceptance_streak = 0;

		if(PRINT) printf("           Best Expectation:   %3.6f  ||  Acceptance Rate: %3.4f (%i/%i)\n",sim_params.best_E,acceptance_rate,proposal_accepted, proposal_count);

		if(poor_acceptance_streak>TEMP_DECAY_LIMIT){
			if(PRINT) printf("NO MC PROGRESS FOR %i TEMP DECAY LIMIT MCDB, TERMINATING\n", TEMP_DECAY_LIMIT);
			break;
		}

		sim_params.temperature=sim_params.temperature*TEMP_EXP_DECAY;
	}
	delete[] k_array, delete[] j_array, delete[] b_array, delete[] k_temp, delete[] j_temp, delete[] b_temp, delete[] sim_params.state, delete[] temp_saved_states;
}





void evolve_mcdb(Simulation_Parameters& sim_params, double *j_array, double *k_array, double *b_array, int random_time_index){
	int i, j;

	for(i=random_time_index; i<sim_params.total_steps; i++){

    if(j_array[i] == 1.0 && k_array[i] == 1.0)	matrix_vector_mult(sim_params.e11,sim_params.state, sim_params.N);
		else if(j_array[i] == 0.0 && k_array[i] == 1.0) matrix_vector_mult(sim_params.e01,sim_params.state, sim_params.N);
		else if(j_array[i] == 1.0 && k_array[i] == 0.0)	matrix_vector_mult(sim_params.e10,sim_params.state, sim_params.N);

		for(j=0; j<2*sim_params.N; j++) sim_params.saved_states[(i+1)*2*sim_params.N +j] = sim_params.state[j];
  }
	if(CHECK) check_norm(sim_params.state, sim_params.N);
}




void calc_initial_temp_mcdb(Simulation_Parameters& sim_params){
	int i, j, random_time_index, start_state_index, count=0;
	double *state_random, *j_temp, *k_temp, *b_temp, sum=0, change, old_E, new_E;

	j_temp           = new double[sim_params.total_steps]();
	k_temp           = new double[sim_params.total_steps]();
	b_temp           = new double[sim_params.total_steps]();
	sim_params.state = new double[2*sim_params.N]();
	state_random     = new double[2*sim_params.N]();

	if(PRINT) printf("\n\n\n\n...Calculating initial temperature based on %i random starting states...\n", RANDOM_STATES);

	for (j=0;j<RANDOM_STATES;j++){
		for (i=0; i<sim_params.N*2;i++) state_random[i] =0.0;
		start_state_index = floor(get_random_double(0,sim_params.N,sim_params.rng));
		state_random[start_state_index*2] = 1;
		std::memcpy(sim_params.state,state_random, 2*sim_params.N*sizeof(double));

		init_arrays_mcdb(sim_params, j_temp, k_temp, b_temp);

		evolve_mcdb(sim_params,j_temp, k_temp, b_temp, 0);
		old_E = cost(sim_params.N, sim_params.state, sim_params.ham_target);

		for (i=0; i<sim_params.total_sweeps*sim_params.sweeps_multiplier;i++){
		      random_time_index = change_array_mcdb(sim_params, j_temp, k_temp, b_temp, i);

			std::memcpy(sim_params.state,state_random, 2*sim_params.N*sizeof(double));//resetting state
			evolve_mcdb(sim_params,j_temp, k_temp, b_temp, 0);
			new_E = cost(sim_params.N, sim_params.state, sim_params.ham_target);

			if (new_E>=old_E) sum += (new_E-old_E), count++;
			old_E=new_E;
		}
	}
	sim_params.temperature = -(sum/(count*log(ACCEPTANCE_PROB)));
	sim_params.initial_temperature = sim_params.temperature;

	delete[] j_temp, delete[] k_temp, delete[] b_temp, delete[] sim_params.state, delete[] state_random;
}

\




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


void scale_best_arrays_mcdb(Simulation_Parameters& sim_params, double* j_best,double* k_best, double* b_best){
	int i;

	for(i=0;i<sim_params.total_steps/2;i++){
		j_best[(sim_params.total_steps-1)-2*i]     = j_best[(sim_params.total_steps/2-1)-i];
		j_best[(sim_params.total_steps-1)-(2*i+1)] = j_best[(sim_params.total_steps/2-1)-i];
		k_best[(sim_params.total_steps-1)-2*i]     = k_best[(sim_params.total_steps/2-1)-i];
		k_best[(sim_params.total_steps-1)-(2*i+1)] = k_best[(sim_params.total_steps/2-1)-i];
		b_best[(sim_params.total_steps-1)-2*i]     = b_best[(sim_params.total_steps/2-1)-i];
		b_best[(sim_params.total_steps-1)-(2*i+1)] = b_best[(sim_params.total_steps/2-1)-i];
	}
}


double change_array_mcdb(Simulation_Parameters& sim_params, double *j_array, double *k_array, double *b_array, int i){
	double *pointer, *pointer2;
	int random_time_index;

	random_time_index = (int)floor(get_random_double(0, sim_params.total_steps, sim_params.rng));


	if(i%2 == 0) pointer = j_array, pointer2 = k_array;
	if(i%2 == 1) pointer = k_array, pointer2 = j_array;
	//if(i%3 == 2) pointer = b_array;

	//rerolling the index if our current index has similar neighbors
	
	if(sim_params.total_steps == MAX_STEPS_MCDB and random_time_index > 0 and random_time_index < sim_params.total_steps-1){
		if(fmod(pointer[random_time_index] + pointer[random_time_index+1] + pointer[random_time_index-1], 3.0) <0.01){
			random_time_index = (int)floor(get_random_double(0, sim_params.total_steps, sim_params.rng));
		}
	}
	
	
	pointer[random_time_index] =  fmod(pointer[random_time_index] + 1.0, 2.0);
	if(pointer[random_time_index]<0.01 and pointer2[random_time_index] <0.01) pointer[random_time_index] = 1;

	return random_time_index;
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


void get_sweeps_mcdb(Simulation_Parameters& sim_params){
	int x = sim_params.total_steps;
	double scalar, a = (double) MAX_STEPS_MCDB, b = (double) MIN_STEPS_MCDB, c = STEPS_CRUNCH_MCDB, offset;

	scalar = (a-b/c)/(a-b);
	offset = b/c - scalar*b;

	sim_params.total_sweeps = (int) ceil(scalar*x + offset) *SWEEPS_MCDB ;

}
