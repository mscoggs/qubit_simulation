#include <gsl/gsl_rng.h>
#include <cstring>
#include <algorithm>    // std::min

#include "check.h"
#include "write_data.h"
#include "hamiltonian.h"
#include "linear_algebra.h"
#include "mcbb.h"
#include "operations.h"
#include "parameters.h"
#include "print.h"




void mcbb_method(Simulation_Parameters& sim_params){
	int i;


	sim_params.init_mcbb_params();

	if(check_commutator(sim_params.N, sim_params.ham_initial, sim_params.ham_target) || sim_params.init_target_dot_squared > INIT_OVERLAP_LIMIT){
		sim_params.tau = 0.0, sim_params.new_distance = 0.0, sim_params.best_mc_result = 0.0;
		if(PRINT) print_mcbb_info(sim_params);
		if(SAVE_DATA) save_mcbb_data_fixed_tau(sim_params);
		sim_params.clear_mcbb_params();
		return;
	}


	while(sim_params.tau<MAX_TAU){

		gsl_rng_set(sim_params.rng, 1);
		calc_initial_temp_mcbb(sim_params);

		for(sim_params.seed=1; sim_params.seed<NUM_SEEDS+1; sim_params.seed++){
			gsl_rng_set(sim_params.rng, sim_params.seed);
			mcbb_simulation(sim_params);

			sim_params.best_mc_result_fixed_tau[sim_params.seed-1] = sim_params.best_mc_result;
			std::memcpy(&sim_params.evolved_state_fixed_tau[(sim_params.seed-1)*2*sim_params.N],sim_params.best_evolved_state, 2*sim_params.N*sizeof(double));
			copy_arrays_mcbb(sim_params,sim_params.j_best_fixed_tau,sim_params.k_best_fixed_tau,sim_params.b_best_fixed_tau,sim_params.j_best,  sim_params.k_best,  sim_params.b_best,  (sim_params.seed-1)*2*NUMBER_OF_BANGS, 0);
		}

		get_best_seed(sim_params);
		if(!update_distances(sim_params)) continue;

		if(PRINT) print_mc_results(sim_params);
		if(SAVE_DATA){
			sim_params.duration = (std::clock() - sim_params.start)/(double) CLOCKS_PER_SEC;
			save_mcbb_data_fixed_tau(sim_params);
		}


		if(sim_params.new_distance < DISTANCE_LIMIT) break;
		else calc_tau(sim_params);
	}
	sim_params.clear_mcbb_params();
}




void mcbb_simulation(Simulation_Parameters& sim_params){
	int i,j, random_time_index, proposal_accepted,proposal_count, poor_acceptance_streak;
	double *j_array, *k_array,*b_array,*j_temp, *k_temp, *b_temp, acceptance_rate, old_mc_result, new_mc_result, change;

	j_array          = new double[2*NUMBER_OF_BANGS]();
	k_array          = new double[2*NUMBER_OF_BANGS]();
	b_array          = new double[2*NUMBER_OF_BANGS]();
	j_temp           = new double[2*NUMBER_OF_BANGS]();
	k_temp           = new double[2*NUMBER_OF_BANGS]();
	b_temp           = new double[2*NUMBER_OF_BANGS]();
	sim_params.state = new double[2*sim_params.N]();


	std::memcpy(sim_params.state,sim_params.init_state, 2*sim_params.N*sizeof(double));
	init_arrays_mcbb(sim_params, sim_params.j_best,sim_params.k_best,sim_params.b_best);
	copy_arrays_mcbb(sim_params, j_array, k_array, b_array,sim_params.j_best,sim_params.k_best,sim_params.b_best,0,0);//a temporary array, used in the undoing of the changes

	evolve_mcbb(sim_params, j_array,k_array,b_array);
	sim_params.best_mc_result = cost(sim_params.target_state, sim_params.N, sim_params.state, sim_params.ham_target);
	old_mc_result = sim_params.best_mc_result;

	if(PRINT) print_mcbb_info(sim_params);

	for (i=0;i<TEMP_DECAY_ITERATIONS;i++){

		proposal_accepted = 0, proposal_count = 0;
		sim_params.temp_iteration = i;

		for (j=0; j<SWEEPS_MCBB*NUMBER_OF_BANGS*sim_params.sweeps_multiplier;j++){

			copy_arrays_mcbb(sim_params, j_temp, k_temp,b_temp,j_array, k_array, b_array,0,0);//a temporary array, used in the undoing of the changes
			change = get_change_mcbb(sim_params);
			random_time_index = (int)floor(get_random_double(0, 2*NUMBER_OF_BANGS, sim_params.rng));
			change_array_mcbb(j_array, k_array, b_array,change,random_time_index,j, sim_params.tau);

			std::memcpy(sim_params.state,sim_params.init_state, 2*sim_params.N*sizeof(double));//resetting state
			evolve_mcbb(sim_params, j_array, k_array, b_array);
			new_mc_result = cost(sim_params.target_state, sim_params.N,sim_params.state, sim_params.ham_target);

			if (new_mc_result<sim_params.best_mc_result) sim_params.best_mc_result=new_mc_result, old_mc_result=new_mc_result,  copy_arrays_mcbb(sim_params, sim_params.j_best, sim_params.k_best, sim_params.b_best,j_array, k_array, b_array, 0, 0), std::memcpy(sim_params.best_evolved_state,sim_params.state, 2*sim_params.N*sizeof(double)), poor_acceptance_streak = 0;
			else if (new_mc_result<=old_mc_result) old_mc_result=new_mc_result;
			else if (get_random_double(0,1,sim_params.rng)<exp(-(new_mc_result-old_mc_result)/(sim_params.temperature))) old_mc_result=new_mc_result, proposal_accepted++, proposal_count++;
			else copy_arrays_mcbb(sim_params, j_array, k_array, b_array,  j_temp, k_temp, b_temp, 0,0), proposal_count++;//undoing the change
		}

		acceptance_rate = (double)proposal_accepted/proposal_count;
		if(acceptance_rate<0.1) poor_acceptance_streak++;
		else poor_acceptance_streak = 0;

		if(PRINT) printf("           Best Expectation:   %3.6f  ||  Acceptance Rate: %3.4f (%i/%i)\n",sim_params.best_mc_result,acceptance_rate,proposal_accepted, proposal_count);

		if(poor_acceptance_streak>TEMP_DECAY_LIMIT){
			if(PRINT) printf("NO MC PROGRESS FOR %i TEMP DECAY LIMIT MCBB, TERMINATING\n", TEMP_DECAY_LIMIT);
			break;
		}

		sim_params.temperature=sim_params.temperature*TEMP_EXP_DECAY;
	}
	delete[] k_array, delete[] j_array, delete[] b_array, delete[] k_temp, delete[] j_temp, delete[] b_temp, delete[] sim_params.state;
}




void evolve_mcbb(Simulation_Parameters& sim_params, double *j_array, double *k_array, double *b_array){
	int *jkb_index, i, j;
	double *hamiltonian,*ham_t_i, *ham_real,*exp_matrix,*e_vals, *v_diag, *jkb_vals, time=0, time2=0, time_step;

	hamiltonian  = new double[2*sim_params.N*sim_params.N]();
	exp_matrix   = new double[2*sim_params.N*sim_params.N]();
	ham_t_i      = new double[2*sim_params.N*sim_params.N]();
	ham_real     = new double[sim_params.N*sim_params.N]();
	v_diag       = new double[sim_params.N*sim_params.N]();
	e_vals       = new double[sim_params.N]();
	jkb_vals     = new double[3]();
	jkb_index    = new int[3]();


	while(time < sim_params.tau){

		time2 = find_next_time(time, sim_params.tau, j_array, k_array, b_array, jkb_index, jkb_vals);
		if(time2 == time) continue; //useless evolution, skipping

		time_step = time2-time;
		time = time2;

		construct_device_hamiltonian_uniform(sim_params, hamiltonian, jkb_vals);
		if(CHECK) check_norm(sim_params.state, sim_params.N);

		if(DIAG){
			for (j=0; j<sim_params.N*sim_params.N; j++) ham_real[j]=0.0,v_diag[j]=0.0, ham_real[j] = hamiltonian[2*j];//converting an all-real-valued complex matrix into just real matrix
			for (j=0; j<sim_params.N; j++) e_vals[j]=0.0;

			diag_hermitian_real_double(sim_params.N, ham_real,v_diag, e_vals);
			exp_diaganolized_real_matrix(hamiltonian, v_diag, e_vals, sim_params.N, time_step);//This function exponentiates e_vals to e^(-i*time_step*e_vals)

			if(CHECK) check_unitary(hamiltonian, sim_params.N);
			matrix_vector_mult(hamiltonian,sim_params.state, sim_params.N);
		}else{
			for (j=0; j<sim_params.N*sim_params.N*2; j++) exp_matrix[j] = 0.0, ham_t_i[j]=0.0;
			for (j=0; j<sim_params.N*sim_params.N; j++) ham_t_i[2*j+1] = (hamiltonian[2*j]*-time_step), ham_t_i[2*j] = (hamiltonian[2*j+1]*time_step);//multiplying by -i*dt for the Pade approximation
			exp_complex_double_matrix_pade(sim_params.N, ham_t_i, exp_matrix);

			if(CHECK) check_unitary(exp_matrix, sim_params.N);
			matrix_vector_mult(exp_matrix, sim_params.state, sim_params.N);
		}
	}
	delete[] hamiltonian, delete[] ham_t_i, delete[] exp_matrix, delete[] e_vals, delete[] v_diag, delete[] ham_real, delete[] jkb_vals, delete[] jkb_index;
}




void calc_initial_temp_mcbb(Simulation_Parameters& sim_params){
	int i, j, random_time_index, init_state_index, count=0;
	double *state_random, *j_temp, *k_temp, *b_temp, sum=0, change, old_mc_result, new_mc_result;

	j_temp           = new double[2*NUMBER_OF_BANGS]();
	k_temp           = new double[2*NUMBER_OF_BANGS]();
	b_temp           = new double[2*NUMBER_OF_BANGS]();
	sim_params.state = new double[2*sim_params.N]();
	state_random     = new double[2*sim_params.N]();


	if(PRINT) printf("\n\n\n\n...Calculating initial temperature based on %i random starting states...\n", RANDOM_STATES);

	for (j=0;j<RANDOM_STATES;j++){
		for (i=0; i<sim_params.N*2;i++) state_random[i] =0.0;
		init_state_index = floor(get_random_double(0,sim_params.N,sim_params.rng));
		state_random[init_state_index*2] = 1;
		std::memcpy(sim_params.state,state_random, 2*sim_params.N*sizeof(double));

		init_arrays_mcbb(sim_params, j_temp, k_temp, b_temp);
		evolve_mcbb(sim_params,j_temp, k_temp, b_temp);
		old_mc_result = cost(sim_params.target_state, sim_params.N, sim_params.state, sim_params.ham_target);

		for (i=0; i<SWEEPS_MCBB*NUMBER_OF_BANGS;i++){
			change = get_change_mcbb(sim_params);
			random_time_index = (int)floor(get_random_double(0, 2*NUMBER_OF_BANGS, sim_params.rng));
			change_array_mcbb(j_temp, k_temp, b_temp,change,random_time_index,j, sim_params.tau);

			std::memcpy(sim_params.state,state_random, 2*sim_params.N*sizeof(double));//resetting state
			evolve_mcbb(sim_params,j_temp, k_temp, b_temp);
			new_mc_result = cost(sim_params.target_state, sim_params.N, sim_params.state, sim_params.ham_target);

			if (new_mc_result>=old_mc_result) sum += (new_mc_result-old_mc_result), count++;
			old_mc_result=new_mc_result;
		}
	}
	sim_params.temperature = -(sum/(count*log(ACCEPTANCE_PROB)));

	delete[] j_temp, delete[] k_temp, delete[] b_temp, delete[] sim_params.state, delete[] state_random;
}





double find_next_time(double time, double tau,double* j_array, double* k_array,double*  b_array,int*  jkb_index, double*jkb){
	bool j = true, k = true, b = true;
	double newtime=tau;


	jkb[0] = fmod(jkb_index[0],2), jkb[1] = fmod(jkb_index[1],2), jkb[2] = fmod(jkb_index[2],2);
	if(jkb_index[0] == NUMBER_OF_BANGS*2) j =false;
	if(jkb_index[1] == NUMBER_OF_BANGS*2) k =false;
	if(jkb_index[2] == NUMBER_OF_BANGS*2) b =false;

	if(!(j || k ||b)) return tau;

	if(j) newtime = j_array[jkb_index[0]];
	if(k  && (k_array[jkb_index[1]] < newtime)) newtime = k_array[jkb_index[1]];
	if(b  && (b_array[jkb_index[2]] < newtime)) newtime = b_array[jkb_index[2]];
	if(newtime == j_array[jkb_index[0]]) jkb_index[0] += 1;
	if(newtime == k_array[jkb_index[1]]) jkb_index[1] += 1;
	if(newtime == b_array[jkb_index[2]]) jkb_index[2] += 1;

	if(jkb_index[0] > 2*NUMBER_OF_BANGS ||  jkb_index[1] > 2*NUMBER_OF_BANGS || jkb_index[2] > 2*NUMBER_OF_BANGS) printf("ERROR: FIND_NEXT_TIME INDEX OUT OF BOUNDS"), exit(0);
	return newtime;
}




void init_arrays_mcbb(Simulation_Parameters& sim_params, double *j_array,double *k_array,double *b_array){
	int i;
	double random_time1, random_time2, random_time3;


	j_array[0] = 0, k_array[0] = 0, b_array[0] = 0;
	j_array[2*NUMBER_OF_BANGS-1] = sim_params.tau, k_array[2*NUMBER_OF_BANGS-1] = sim_params.tau, b_array[2*NUMBER_OF_BANGS-1] = sim_params.tau;
	for (i=1; i<NUMBER_OF_BANGS*2-1;i++){
		random_time1 = get_random_double(0,TAU_INIT,sim_params.rng);
	       	random_time2 = get_random_double(0,TAU_INIT,sim_params.rng);
	       	random_time3 = get_random_double(0,TAU_INIT,sim_params.rng);
		j_array[i] = random_time1, k_array[i] = random_time2, b_array[i] = random_time3;
	}
	std::sort(j_array, j_array + NUMBER_OF_BANGS*2), std::sort(k_array, k_array + NUMBER_OF_BANGS*2), std::sort(b_array, b_array + NUMBER_OF_BANGS*2);
}




void copy_arrays_mcbb(Simulation_Parameters& sim_params, double* j_to,  double* k_to, double* b_to, double *j_from, double *k_from, double* b_from, int destination_index, int source_index){


	std::memcpy(&j_to[destination_index], &j_from[source_index], 2*NUMBER_OF_BANGS*sizeof(double));
	std::memcpy(&k_to[destination_index], &k_from[source_index], 2*NUMBER_OF_BANGS*sizeof(double));
	std::memcpy(&b_to[destination_index], &b_from[source_index], 2*NUMBER_OF_BANGS*sizeof(double));
}




void change_array_mcbb(double *j_array, double *k_array, double *b_array, double change, int random_time_index, int i, double tau){
	double upper_bound , lower_bound, *pointer;


	if(i%2 == 0) pointer = j_array;
	if(i%2 == 1) pointer = k_array;
	//if(i%3 == 2) pointer = b_array;

	if(change > 0){
		if(random_time_index == 2*NUMBER_OF_BANGS - 1) upper_bound = tau;
		else upper_bound = *(pointer + random_time_index +1);

		if( *(pointer+random_time_index) + change < upper_bound) *(pointer+random_time_index) += change;
		else *(pointer+random_time_index) = upper_bound;
	}else{
		if(random_time_index == 0) lower_bound = 0;
		else lower_bound = *(pointer+random_time_index-1);

		if( *(pointer+random_time_index) + change > lower_bound) *(pointer+random_time_index) += change;
		else *(pointer+random_time_index) = lower_bound;
	}
}




double get_change_mcbb(Simulation_Parameters& sim_params){
	int sign;
	double max_change, abs_change, temp_scalar;

	temp_scalar = sim_params.temp_iteration/(double)TEMP_DECAY_ITERATIONS;

	sign       = pow(-1,(int)floor(get_random_double(0,10,sim_params.rng)));
	max_change = sim_params.tau*(MAX_CHANGE_FRACTION_MCBB*(1-temp_scalar) + MIN_CHANGE_FRACTION_MCBB*temp_scalar);
	abs_change = get_random_double(0,max_change, sim_params.rng);

	return sign*abs_change;
}




void scale_arrays_mcbb(double *j_array, double *k_array, double *b_array, double scalar){
	int i;


	for(i=0;i<NUMBER_OF_BANGS*2;i++) j_array[i] = j_array[i]*scalar ,k_array[i] = k_array[i]*scalar, b_array[i] = b_array[i]*scalar;
}








void binary_search_mcbb(Simulation_Parameters& sim_params){
	// int i;
	// double tau_max = sim_params.tau, tau_min, distance_temp;
	//
	//
	// tau_min = sim_params.tau/TAU_SCALAR_TINY;
	// sim_params.tau = (tau_max + tau_min) / 2.0;
	// sim_params.new_distance = sim_params.old_distance;
	//
	// if(PRINT) printf("\nUsing binary search method to look for optimal ground state...");
	//
	// while((tau_max - tau_min) >BINARY_SEARCH_TAU_LIMIT_MCBB){
	//
	// 	gsl_rng_set(sim_params.rng, 1);
	// 	calc_initial_temp_mcbb(sim_params);
	//
	// 	for(sim_params.seed=1; sim_params.seed<NUM_SEEDS+1; sim_params.seed++){
	// 		gsl_rng_set(sim_params.rng, sim_params.seed);
	// 		mcbb_simulation(sim_params);
	//
	// 		sim_params.best_mc_result_fixed_tau[sim_params.seed-1] = sim_params.best_mc_result;
	// 		copy_arrays_mcbb(sim_params,sim_params.j_best_fixed_tau,sim_params.k_best_fixed_tau,sim_params.b_best_fixed_tau,sim_params.j_best,  sim_params.k_best,  sim_params.b_best,  (sim_params.seed-1)*2*NUMBER_OF_BANGS, 0);
	// 	}
	// 	for(sim_params.seed=1; sim_params.seed<NUM_SEEDS+1; sim_params.seed++) if(sim_params.best_mc_result_fixed_tau[sim_params.seed-1] < sim_params.best_mc_result) sim_params.best_mc_result = sim_params.best_mc_result_fixed_tau[sim_params.seed-1];
	//
	// 	distance_temp	= sim_params.new_distance;
	// 	sim_params.new_distance = calc_distance(sim_params);
	//
	// 	if(PRINT) print_mc_results(sim_params);
	//
	// 	if(sim_params.new_distance < DISTANCE_LIMIT){
	// 		if(PRINT) printf("\nIn binary_search_mcbb....\nStepping backward....\n");
	// 		tau_max = sim_params.tau;
	// 		sim_params.tau = (tau_max+tau_min)/2.0;
	//
	// 	}else{
	// 		if(PRINT) printf("\nIn binary_search_mcbb....\nStepping forward....\n");
	// 		sim_params.old_distance = sim_params.new_distance;
	// 		if(SAVE_DATA) save_mcbb_data_fixed_tau(sim_params);
	//
	// 		tau_min = sim_params.tau;
	// 		sim_params.tau = (tau_min+tau_max)/2;
	// 	}
	// }
	// if(PRINT) printf("\nUpper and lower tau are close enough, exiting binary_search_mcbb\n");
}



void evolve_fixed_protocol(Simulation_Parameters& sim_params){
	// int i;
	// double new_mc_result, j_fraction_start = 0.0, j_fraction_end = 1.0, k_fraction_start = 0.52, k_fraction_end = 1.0;
	//
	//
	// sim_params.init_mcbb_params();
	// sim_params.state                = new double[2*sim_params.N]();
	//
	// std::memcpy(sim_params.state,sim_params.init_state, 2*sim_params.N*sizeof(double));
	//

	//
	// while(sim_params.tau<1.0){
	// 	for(i=0;i<2*NUMBER_OF_BANGS;i++) sim_params.j_best[i] = sim_params.tau, sim_params.k_best[i] = sim_params.tau;
	// 	sim_params.j_best[0] =	j_fraction_start*sim_params.tau;
	// 	sim_params.k_best[0] =	k_fraction_start*sim_params.tau;
	// 	sim_params.j_best[1] =	j_fraction_end*sim_params.tau;
	// 	sim_params.k_best[1] =	k_fraction_end*sim_params.tau;
	//
	// 	std::memcpy(sim_params.state,sim_params.init_state, 2*sim_params.N*sizeof(double));
	// 	evolve_mcbb(sim_params, sim_params.j_best, sim_params.k_best, sim_params.b_best);
	// 	new_mc_result = cost(sim_params.N,sim_params.state, sim_params.ham_target);
	// 	printf("tau: %f, E: %f\n", sim_params.tau, new_mc_result);
	// 	sim_params.tau = sim_params.tau*1.05;
	// }
	// sim_params.clear_mcbb_params();
}
