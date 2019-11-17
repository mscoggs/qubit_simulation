#include <gsl/gsl_rng.h>
#include <cstring>
#include <algorithm>    // std::min

#include "check.h"
#include "export_data.h"
#include "hamiltonian.h"
#include "linear_algebra.h"
#include "mcbb.h"
#include "operations.h"
#include "parameters.h"
#include "print.h"



void mcbb_method(Simulation_Parameters& sim_params){
	int i;
	double distance_temp;
	sim_params.best_arrays_size = ceil(2*NUMBER_OF_BANGS*(MAX_TAU_STEPS_MC+MAX_BS_STEPS_MC));
	sim_params.tau_array_size   = ceil(MAX_TAU_STEPS_MC+MAX_BS_STEPS_MC);
	sim_params.j_best 					= new double[sim_params.best_arrays_size]();
	sim_params.k_best 					= new double[sim_params.best_arrays_size]();
	sim_params.b_best 					= new double[sim_params.best_arrays_size]();
	sim_params.tau_array 				= new double[sim_params.tau_array_size]();
	sim_params.best_E_array		  = new double[sim_params.tau_array_size]();

	for(sim_params.seed=1; sim_params.seed<SEED_TOTAL+1 ;sim_params.seed++){

		gsl_rng_set(sim_params.rng, sim_params.seed);
		for (i=0;i<sim_params.best_arrays_size;i++) sim_params.j_best[i] = 0,sim_params.k_best[i] = 0,sim_params.b_best[i] = 0;
		for (i=0;i<sim_params.tau_array_size;i++) sim_params.best_E_array[i] = 0, sim_params.tau_array[i] = 0;


		sim_params.old_distance = 1;
		sim_params.new_distance =1;

		sim_params.tau = TAU_INIT_MCBB;
		sim_params.tau_old = TAU_INIT_MCBB;
		sim_params.tau_array[0] = sim_params.tau;

		sim_params.source_index = 0;
		sim_params.destination_index = 2*NUMBER_OF_BANGS;
		sim_params.index = 1;

		init_arrays_mcbb(sim_params, sim_params.j_best, sim_params.k_best, sim_params.b_best);


		while(sim_params.tau<MAX_TAU_MCBB){

			calc_initial_temp_mcbb(sim_params);
			mcbb_simulation(sim_params);
			distance_temp = calc_distance(sim_params.best_E_array[0], sim_params.ground_E,  sim_params.best_E);
			if(PRINT) print_mc_results(sim_params);

			if(distance_temp < DIFFERENCE_LIMIT_MCBB){
				binary_search_mcbb(sim_params);
				break;

			}else{

				sim_params.old_distance = sim_params.new_distance;
				sim_params.new_distance = distance_temp;

				sim_params.tau_array[sim_params.index] = sim_params.tau;
				sim_params.best_E_array[sim_params.index] = sim_params.best_E;

				sim_params.tau_old = sim_params.tau;
				calc_tau_mcbb(sim_params);

				sim_params.source_index = sim_params.destination_index;
				sim_params.destination_index = sim_params.destination_index + 2*NUMBER_OF_BANGS;

				sim_params.index ++;
			}
		}
		sim_params.tau = sim_params.tau_old;
		if(MCBB_DATA) save_mcbb_data(sim_params);
	}
	delete[] sim_params.j_best, delete[] sim_params.k_best, delete[] sim_params.b_best, delete[] sim_params.tau_array, delete[] sim_params.best_E_array;
}



void mcbb_simulation(Simulation_Parameters& sim_params){
	int i,j, random_time_index, proposal_accepted,proposal_count, poor_acceptance_streak=0;
	double *j_array, *k_array,*b_array,*j_temp, *k_temp, *b_temp, acceptance_rate, old_E, new_E, change;

	j_array = new double[2*NUMBER_OF_BANGS]();
	k_array = new double[2*NUMBER_OF_BANGS]();
	b_array = new double[2*NUMBER_OF_BANGS]();
	j_temp = new double[2*NUMBER_OF_BANGS]();
	k_temp = new double[2*NUMBER_OF_BANGS]();
	b_temp = new double[2*NUMBER_OF_BANGS]();
	sim_params.state =  new double[2*sim_params.N]();
	std::memcpy(sim_params.state,sim_params.start_state, 2*sim_params.N*sizeof(double));


	copy_arrays_mcbb(sim_params, j_array,  k_array,  b_array,sim_params.j_best,  sim_params.k_best,  sim_params.b_best,  0, sim_params.source_index);
	scale_arrays_mcbb(j_array,k_array,b_array,sim_params.tau/sim_params.tau_old);
	copy_arrays_mcbb(sim_params, sim_params.j_best,  sim_params.k_best,  sim_params.b_best,j_array,k_array,b_array,   sim_params.destination_index, 0);//This is necessary to cover the case that this first scaled array is the best protocol. If this weren't here, this array would never get appended in the iteration below

	evolve_mcbb(sim_params, j_array,k_array,b_array);
	sim_params.best_E = cost(sim_params.N, sim_params.state, sim_params.ham_target);
	old_E = sim_params.best_E;

	if(sim_params.index == 1) sim_params.best_E_array[0] = sim_params.best_E;
	if(PRINT) print_mcbb_info(sim_params), printf("Pre-Monte_Carlo Expectation:   %f\n", sim_params.best_E);

	for (i=0;i<TEMP_DECAY_ITERATIONS_MCBB;i++){
		proposal_accepted = 0, proposal_count = 0;

		for (j=0; j<SWEEPS_MCBB*NUMBER_OF_BANGS;j++){

			copy_arrays_mcbb(sim_params, j_temp, k_temp,b_temp,j_array, k_array, b_array,0,0);//a temporary array, used in the undoing of the changes

			change = get_change_mcbb(sim_params);
			random_time_index = (int)floor(get_random_double(0, 2*NUMBER_OF_BANGS, sim_params.rng));

			change_array_mcbb(j_array, k_array, b_array,change,random_time_index,j, sim_params.tau);
			std::memcpy(sim_params.state,sim_params.start_state, 2*sim_params.N*sizeof(double));//resetting state

			evolve_mcbb(sim_params, j_array, k_array, b_array);
			new_E = cost(sim_params.N,sim_params.state, sim_params.ham_target);

			if (new_E<sim_params.best_E) sim_params.best_E=new_E, old_E=new_E,  copy_arrays_mcbb(sim_params, sim_params.j_best, sim_params.k_best, sim_params.b_best,j_array, k_array, b_array, sim_params.destination_index, 0), poor_acceptance_streak = 0;
			else if (new_E<=old_E) old_E=new_E;
			else if (get_random_double(0,1,sim_params.rng)<exp(-(new_E-old_E)/(sim_params.temperature))) old_E=new_E, proposal_accepted++, proposal_count++;
			else copy_arrays_mcbb(sim_params, j_array, k_array, b_array,  j_temp, k_temp, b_temp, 0,0), proposal_count++;//undoing the change
		}

		acceptance_rate = (double)proposal_accepted/proposal_count;
		if(acceptance_rate<0.1) poor_acceptance_streak++;
		else poor_acceptance_streak = 0;

		if(PRINT) printf("           Best Expectation:   %3.6f  ||  Acceptance Rate: %3.4f (%i/%i)\n",sim_params.best_E,acceptance_rate,proposal_accepted, proposal_count);

		if(poor_acceptance_streak>TEMP_DECAY_LIMIT_MCBB){
			printf("NO MC PROGRESS FOR %i TEMP DECAY LIMIT MCBB, TERMINATING\n", TEMP_DECAY_LIMIT_MCBB);
			break;
		}
		sim_params.temperature=sim_params.temperature*TEMP_EXP_DECAY_MC;
	}
	delete[] k_array, delete[] j_array, delete[] b_array, delete[] k_temp, delete[] j_temp, delete[] b_temp, delete[] sim_params.state;
}



void evolve_mcbb(Simulation_Parameters& sim_params, double *j_array, double *k_array, double *b_array){
	int *jkb_index, i, j;
	double *hamiltonian,*ham_t_i, *ham_real,*exp_matrix,*e_vals, *v_diag, *jkb_vals, time=0, time2=0, time_step;
	hamiltonian    		 = new double[2*sim_params.N*sim_params.N]();
	exp_matrix = new double[2*sim_params.N*sim_params.N]();
	ham_t_i    = new double[2*sim_params.N*sim_params.N]();
	ham_real   = new double[sim_params.N*sim_params.N]();
	v_diag     = new double[sim_params.N*sim_params.N]();
	e_vals     = new double[sim_params.N]();
	jkb_index  = new int[3]();
	jkb_vals   = new double[3]();

	while(time < sim_params.tau){
		time2 = find_next_time(time, sim_params.tau, j_array, k_array, b_array, jkb_index, jkb_vals);
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
	if(PRINT) printf("\n\n\n\n...Calculating initial temperature based on %i random starting states...\n", RANDOM_STATES_MCBB);

	int i, j, random_time_index, start_state_index=0, count=0;
	double *state_random, *j_temp, *k_temp, *b_temp, sum=0, change, old_E, new_E;

	j_temp = new double[2*NUMBER_OF_BANGS]();
	k_temp = new double[2*NUMBER_OF_BANGS]();
	b_temp = new double[2*NUMBER_OF_BANGS]();
	sim_params.state =  new double[2*sim_params.N]();
	state_random = new double[2*sim_params.N]();

	for (j=0;j<RANDOM_STATES_MCBB;j++){

		for (i=0; i<sim_params.N*2;i++) state_random[i] =0.0;
		start_state_index = floor(get_random_double(0,sim_params.N,sim_params.rng));
		state_random[start_state_index*2] = 1;
		std::memcpy(sim_params.state,state_random, 2*sim_params.N*sizeof(double));

		init_arrays_mcbb(sim_params, j_temp, k_temp, b_temp);
		evolve_mcbb(sim_params,j_temp, k_temp, b_temp);
		old_E = cost(sim_params.N, sim_params.state, sim_params.ham_target);

		for (i=0; i<SWEEPS_MCBB*NUMBER_OF_BANGS;i++){
			std::memcpy(sim_params.state,state_random, 2*sim_params.N*sizeof(double));//resetting state
			change = get_change_mcbb(sim_params);
			random_time_index = (int)floor(get_random_double(0, 2*NUMBER_OF_BANGS, sim_params.rng));

			change_array_mcbb(j_temp, k_temp, b_temp,change,random_time_index,j, sim_params.tau);
			evolve_mcbb(sim_params,j_temp, k_temp, b_temp);
			new_E = cost(sim_params.N, sim_params.state, sim_params.ham_target);

			if (new_E>=old_E) sum += (new_E-old_E), count++;
			old_E=new_E;
		}
	}

	delete[] j_temp, delete[] k_temp, delete[] b_temp, delete[] sim_params.state, delete[] state_random;
	sim_params.temperature = -(sum/(count*log(ACCEPTANCE_PROB_MC)));
}



void binary_search_mcbb(Simulation_Parameters& sim_params){
	printf("\nUsing binary search method to look for optimal ground state...");

	double tau_max = sim_params.tau, tau_min = sim_params.tau_array[sim_params.index-1], distance_temp;
	sim_params.tau = (tau_max + tau_min) / 2.0;

	while((tau_max - tau_min) >BINARY_SEARCH_TAU_LIMIT_MCBB){

		calc_initial_temp_mcbb(sim_params);
		mcbb_simulation(sim_params);
		distance_temp = calc_distance(sim_params.best_E_array[0], sim_params.ground_E,  sim_params.best_E);

		if(PRINT) print_mc_results(sim_params);

		if(distance_temp < DIFFERENCE_LIMIT_MCBB){//still in the ground state, backtrack tau
			printf("\nIn binary_search_mcbf....\nStepping backward....\n");

			tau_max = sim_params.tau;
			sim_params.tau_old = sim_params.tau;
			sim_params.tau = (tau_max+tau_min)/2.0;

		}else{
			printf("\nIn binary_search_mcbf....\nStepping forward....\n");

			sim_params.old_distance = sim_params.new_distance;
			sim_params.new_distance = distance_temp;

			sim_params.tau_array[sim_params.index] = sim_params.tau;
			sim_params.best_E_array[sim_params.index] = sim_params.best_E;

			tau_min = sim_params.tau;
			sim_params.tau_old = sim_params.tau;
			sim_params.tau = (tau_min+tau_max)/2;

			sim_params.source_index = sim_params.destination_index;
			sim_params.destination_index = sim_params.destination_index + 2*NUMBER_OF_BANGS;

			sim_params.index ++;
		}
	}
	printf("\nUpper and lower tau are close enough, exiting binary_search_mcbb\n");
}



void calc_tau_mcbb(Simulation_Parameters& sim_params){
	double difference = abs(sim_params.old_distance - sim_params.new_distance), tau_scalar;

	if((difference < .1) and (sim_params.new_distance > 0.2 )) tau_scalar = TAU_SCALAR_MC * 2 ;
	else tau_scalar = TAU_SCALAR_MCBB;
	sim_params.tau = sim_params.tau*tau_scalar;
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

	if(jkb_index[0] > 2*NUMBER_OF_BANGS ||  jkb_index[0] > 2*NUMBER_OF_BANGS || jkb_index[0] > 2*NUMBER_OF_BANGS) printf("ERROR: FIND_NEXT_TIME INDEX OUT OF BOUNDS"), exit(0);

	return newtime;
}



void init_arrays_mcbb(Simulation_Parameters& sim_params, double *j_array,double *k_array,double *b_array){
	int i;
	double random_time_index1, random_time_index2, random_time_index3,  tau = TAU_INIT_MCBB;
	j_array[0] = 0, k_array[0] = 0, b_array[0] = 0;
	j_array[2*NUMBER_OF_BANGS-1] = tau, k_array[2*NUMBER_OF_BANGS-1] = tau, b_array[2*NUMBER_OF_BANGS-1] = tau;
	for (i=1; i<NUMBER_OF_BANGS*2-1;i++){
		random_time_index1 = get_random_double(0,TAU_INIT_MCBB,sim_params.rng), random_time_index2 = get_random_double(0,TAU_INIT_MCBB,sim_params.rng), random_time_index3 = get_random_double(0,TAU_INIT_MCBB,sim_params.rng);
		j_array[i] = random_time_index1, k_array[i] = random_time_index2, b_array[i] = random_time_index3;
	}
	std::sort(j_array, j_array + NUMBER_OF_BANGS*2), std::sort(k_array, k_array + NUMBER_OF_BANGS*2), std::sort(b_array, b_array + NUMBER_OF_BANGS*2);
}



void copy_arrays_mcbb(Simulation_Parameters& sim_params, double* j_to,  double* k_to, double* b_to, double *j_from, double *k_from, double* b_from, int destination_index, int source_index){
	std::memcpy(&j_to[destination_index], &j_from[source_index], 2*NUMBER_OF_BANGS*sizeof(double));
	std::memcpy(&k_to[destination_index], &k_from[source_index], 2*NUMBER_OF_BANGS*sizeof(double));
	std::memcpy(&b_to[destination_index], &b_from[source_index], 2*NUMBER_OF_BANGS*sizeof(double));
}



void change_array_mcbb(double *j_array, double *k_array, double *b_array, double change, int random_time_index, int i, double tau){
	double upper_bound, lower_bound, *pointer;
	if(i%3 == 0) pointer = j_array;
	if(i%3 == 1) pointer = k_array;
	if(i%3 == 2) pointer = b_array;

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
	int sign = pow(-1,(int)floor(get_random_double(0,10,sim_params.rng)));
	double max_change = sim_params.tau*CHANGE_FRACTION_MCBB;
  double abs_change = get_random_double(0,max_change, sim_params.rng);
	return sign*abs_change;
}



void scale_arrays_mcbb(double *j_array, double *k_array, double *b_array, double scalar){
	int i;
	for(i=0;i<NUMBER_OF_BANGS*2;i++) j_array[i] = j_array[i]*scalar ,k_array[i] = k_array[i]*scalar, b_array[i] = b_array[i]*scalar;
}
