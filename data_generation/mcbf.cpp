#include <gsl/gsl_rng.h>
#include <cstring>
#include <algorithm>    // std::min

#include "check.h"
#include "export_data.h"
#include "hamiltonian.h"
#include "linear_algebra.h"
#include "mcbf.h"
#include "operations.h"
#include "parameters.h"
#include "print.h"



void mcbf_method(Simulation_Parameters& sim_params){
	int i;
	sim_params.init_mcbf_params();

	if(check_commutator(sim_params.N, sim_params.ham_initial, sim_params.ham_target) || sim_params.initial_E -sim_params.ground_E < 0.001){
		sim_params.tau = 0.0, sim_params.new_distance = 0.0, sim_params.best_E = 0.0;
		if(SAVE_DATA) save_mcbf_data_fixed_tau(sim_params);
		sim_params.clear_mcbf_params();
		return;
	}

	while(sim_params.tau<MAX_TAU){

		gsl_rng_set(sim_params.rng, 1);
		calc_initial_temp_mcbf(sim_params);

		for(sim_params.seed=1; sim_params.seed<NUM_SEEDS+1; sim_params.seed++){
			gsl_rng_set(sim_params.rng, sim_params.seed);
			mcbf_simulation(sim_params);

			sim_params.E_array_fixed_tau[sim_params.seed-1] = sim_params.best_E;
			std::memcpy(&sim_params.evolved_state_fixed_tau[(sim_params.seed-1)*2*sim_params.N],sim_params.best_evolved_state, 2*sim_params.N*sizeof(double));
			copy_arrays_mcbf(sim_params,sim_params.j_best_fixed_tau,sim_params.k_best_fixed_tau,sim_params.b_best_fixed_tau,sim_params.j_best,  sim_params.k_best,  sim_params.b_best,  (sim_params.seed-1)*NUMBER_OF_SITES*sim_params.total_steps, 0);
		}
		get_best_seed(sim_params);

		if(!update_distances(sim_params)) continue;

		if(PRINT) print_mc_results(sim_params);
		if(SAVE_DATA){
			sim_params.duration = (std::clock() - sim_params.start)/(double) CLOCKS_PER_SEC;
		       	save_mcbf_data_fixed_tau(sim_params);
		}

		if(sim_params.new_distance < DISTANCE_LIMIT) break;
		else{
			calc_tau(sim_params);
			calc_total_steps_mcbf(sim_params);
			sim_params.time_step = sim_params.tau/((double) sim_params.total_steps);
		}
	}
	sim_params.clear_mcbf_params();
}


void mcbf_simulation(Simulation_Parameters& sim_params){
	int i,j, random_row, random_col, proposal_accepted,proposal_count, poor_acceptance_streak;
	double *j_array, *k_array,*b_array,*j_temp, *k_temp, *b_temp, acceptance_rate, old_E, new_E, change;

	j_array          = new double[2*NUMBER_OF_SITES*sim_params.total_steps]();
	k_array          = new double[2*NUMBER_OF_SITES*sim_params.total_steps]();
	b_array          = new double[NUMBER_OF_SITES*sim_params.total_steps]();
	j_temp           = new double[2*NUMBER_OF_SITES*sim_params.total_steps]();
	k_temp           = new double[2*NUMBER_OF_SITES*sim_params.total_steps]();
	b_temp           = new double[NUMBER_OF_SITES*sim_params.total_steps]();
	sim_params.state = new double[2*sim_params.N]();


	std::memcpy(sim_params.state,sim_params.init_state, 2*sim_params.N*sizeof(double));
	init_arrays_mcbf(sim_params, sim_params.j_best,sim_params.k_best,sim_params.b_best);
	copy_arrays_mcbf(sim_params, j_array, k_array, b_array,sim_params.j_best,sim_params.k_best,sim_params.b_best,0,0);//a temporary array, used in the undoing of the changes

	evolve_mcbf(sim_params, j_array,k_array,b_array);
	sim_params.best_E = cost(sim_params.N, sim_params.state, sim_params.ham_target);
	old_E = sim_params.best_E;

	if(PRINT) print_mcbf_info(sim_params);

	for (i=0;i<TEMP_DECAY_ITERATIONS;i++){

		proposal_accepted = 0,proposal_count = 0;

		for (j=0; j<SWEEPS_MC*sim_params.total_steps*sim_params.sweeps_multiplier;j++){

			copy_arrays_mcbf(sim_params, j_temp, k_temp,b_temp,j_array, k_array, b_array,0,0);//a temporary array, used in the undoing of the changes

			change = get_change_mcbf(sim_params);
			random_row = floor(get_random_double(0,sim_params.total_steps,sim_params.rng));
			random_col = floor(get_random_double(0,NUMBER_OF_SITES*2,sim_params.rng));
			change_array_mcbf(j_array,k_array,b_array,random_row,random_col,change,j, sim_params.total_steps);

			std::memcpy(sim_params.state,sim_params.init_state, 2*sim_params.N*sizeof(double));//resetting state
			evolve_mcbf(sim_params, j_array,k_array,b_array);
			new_E = cost(sim_params.N,sim_params.state, sim_params.ham_target);

			if (new_E<sim_params.best_E) sim_params.best_E=new_E, old_E=new_E,  copy_arrays_mcbf(sim_params, sim_params.j_best, sim_params.k_best, sim_params.b_best,j_array, k_array, b_array, 0, 0), std::memcpy(sim_params.best_evolved_state,sim_params.state, 2*sim_params.N*sizeof(double)), poor_acceptance_streak = 0;
			else if (new_E<=old_E) old_E=new_E;
			else if (get_random_double(0,1,sim_params.rng)<exp(-(new_E-old_E)/(sim_params.temperature))) old_E=new_E, proposal_accepted++, proposal_count++;
			else copy_arrays_mcbf(sim_params, j_array, k_array, b_array,  j_temp, k_temp, b_temp, 0,0), proposal_count++;//undoing the change
		}

		acceptance_rate = (double)proposal_accepted/proposal_count;
		if(acceptance_rate<0.1) poor_acceptance_streak++;
		else poor_acceptance_streak = 0;

		if(PRINT) printf("           Best Expectation:   %3.6f  ||  Acceptance Rate: %3.4f (%i/%i)\n",sim_params.best_E,acceptance_rate,proposal_accepted, proposal_count);

		if(poor_acceptance_streak>TEMP_DECAY_LIMIT){
			if(PRINT) printf("NO MC PROGRESS FOR %i TEMP DECAY ITERATIONS, TERMINATING\n", TEMP_DECAY_LIMIT);
			break;
		}

		sim_params.temperature=sim_params.temperature*TEMP_EXP_DECAY;
	}
	delete[] j_array, delete[] k_array, delete[] b_array, delete[] j_temp, delete[] k_temp, delete[] b_temp, delete[] sim_params.state;
}



void evolve_mcbf(Simulation_Parameters& sim_params, double *j_array, double *k_array, double *b_array){
	int i,j;
	double *hamiltonian,*ham_t_i, *ham_real,*exp_matrix,*e_vals, *v_diag;

	hamiltonian = new double[2*sim_params.N*sim_params.N]();
	exp_matrix  = new double[2*sim_params.N*sim_params.N]();
	ham_t_i     = new double[2*sim_params.N*sim_params.N]();
	ham_real    = new double[sim_params.N*sim_params.N]();
	v_diag      = new double[sim_params.N*sim_params.N]();
	e_vals      = new double[sim_params.N]();

	for (i=0; i<sim_params.total_steps;i++){
		construct_device_hamiltonian(sim_params, hamiltonian, j_array, k_array, b_array,i);

		if(DIAG){
			for (j=0; j<sim_params.N*sim_params.N; j++) ham_real[j]=0.0,v_diag[j]=0.0,ham_real[j] = hamiltonian[2*j];  //converting an all-real-valued complex matrix into just real matrix
			for (j=0; j<sim_params.N; j++) e_vals[j]=0.0;

			diag_hermitian_real_double(sim_params.N, ham_real,v_diag, e_vals);
			exp_diaganolized_real_matrix(hamiltonian, v_diag, e_vals, sim_params.N, sim_params.time_step);//This function exponentiates e_vals to e^(-i*time_step*e_vals)

			if(CHECK) check_unitary(hamiltonian, sim_params.N);
			matrix_vector_mult(hamiltonian,sim_params.state, sim_params.N);
		}
		else{
			for (j=0; j<sim_params.N*sim_params.N*2; j++) exp_matrix[j] = 0.0, ham_t_i[j]=0.0;
			for (j=0; j<sim_params.N*sim_params.N; j++) ham_t_i[2*j+1] = (hamiltonian[2*j]*-sim_params.time_step), ham_t_i[2*j] = (hamiltonian[2*j+1]*sim_params.time_step);  //multiplying by -i*dt for the Pade approximation
			exp_complex_double_matrix_pade(sim_params.N, ham_t_i, exp_matrix);

			if(CHECK) check_unitary(exp_matrix, sim_params.N);
			matrix_vector_mult(exp_matrix, sim_params.state, sim_params.N);
		}
		if(CHECK) check_norm(sim_params.state, sim_params.N);
	}
	delete[] hamiltonian, delete[] ham_t_i, delete[] exp_matrix, delete[] e_vals, delete[] v_diag, delete[] ham_real;
}




void calc_initial_temp_mcbf(Simulation_Parameters& sim_params){
	int i, j, random_row, random_col, init_state_index, count=0;
	double *state_random,*j_temp, *k_temp, *b_temp, sum=0, change, old_E, new_E;

	j_temp           = new double[2*NUMBER_OF_SITES*sim_params.total_steps]();
	k_temp           = new double[2*NUMBER_OF_SITES*sim_params.total_steps]();
	b_temp           = new double[NUMBER_OF_SITES*sim_params.total_steps]();
	sim_params.state = new double[2*sim_params.N]();
	state_random     = new double[2*sim_params.N]();


	if(PRINT) printf("\n\n\n\n...Calculating initial temperature based on %i random starting states...\n", RANDOM_STATES);

	for (i=0;i<RANDOM_STATES;i++){
		for (j=0; j<sim_params.N*2;j++) state_random[j] =0.0;
		init_state_index = (int) floor(get_random_double(0,sim_params.N,sim_params.rng));
		state_random[init_state_index*2] = 1;
		std::memcpy(sim_params.state,state_random, 2*sim_params.N*sizeof(double));

		init_arrays_mcbf(sim_params, j_temp, k_temp, b_temp);
		evolve_mcbf(sim_params, j_temp,k_temp,b_temp);
		old_E = cost(sim_params.N, sim_params.state, sim_params.ham_target);

		for (j=0; j<SWEEPS_MC*sim_params.total_steps;j++){
			change = get_change_mcbf(sim_params);
			random_row = floor(get_random_double(0,sim_params.total_steps,sim_params.rng));
			random_col = floor(get_random_double(0,NUMBER_OF_SITES*2,sim_params.rng));
			change_array_mcbf(j_temp,k_temp,b_temp,random_row,random_col,change,j, sim_params.total_steps);

			std::memcpy(sim_params.state,state_random, 2*sim_params.N*sizeof(double));//resetting state
			evolve_mcbf(sim_params, j_temp,k_temp,b_temp);
			new_E = cost(sim_params.N, sim_params.state, sim_params.ham_target);
			if (new_E>old_E) sum += (new_E-old_E), count++;
			old_E=new_E;
		}
	}
	sim_params.temperature = -(sum/(count*log(ACCEPTANCE_PROB)));

	delete[] j_temp, delete[] k_temp, delete[] b_temp, delete[] state_random, delete[] sim_params.state;
}






void calc_total_steps_mcbf(Simulation_Parameters& sim_params){


	sim_params.total_steps = ceil((1-sim_params.new_distance)*MAX_EVOLVE_STEPS_MC + sim_params.new_distance*TOTAL_STEPS_INIT_MC);
}



void init_arrays_mcbf(Simulation_Parameters& sim_params, double *j_array, double *k_array, double *b_array){
	int i,j;
	double random_val1, random_val2, random_val3;


	if(UNIFORM_SITES){
		for (i=0; i<sim_params.total_steps;i++){
			random_val1 = get_random_double(MIN_PARAM,MAX_PARAM,sim_params.rng);
		       	random_val2 = get_random_double(MIN_PARAM,MAX_PARAM,sim_params.rng);
		       	random_val3 = get_random_double(MIN_PARAM,MAX_PARAM,sim_params.rng);
			for (j=0; j<NUMBER_OF_SITES*2; j++){
			       	j_array[i*NUMBER_OF_SITES*2+j] = random_val2;
			       	k_array[i*NUMBER_OF_SITES*2+j] = random_val3;
			}
			for (j=0; j<NUMBER_OF_SITES; j++) b_array[i*NUMBER_OF_SITES+j]= random_val1;
		}
	}else{
		for (i=0; i<sim_params.total_steps;i++){
			for (j=0; j<NUMBER_OF_SITES*2; j++){
			       	j_array[i*NUMBER_OF_SITES*2+j] = get_random_double(MIN_PARAM,MAX_PARAM,sim_params.rng);
			       	k_array[i*NUMBER_OF_SITES*2+j] = get_random_double(MIN_PARAM,MAX_PARAM,sim_params.rng);
			}
			for (j=0; j<NUMBER_OF_SITES; j++) b_array[i*NUMBER_OF_SITES+j]= get_random_double(MIN_PARAM,MAX_PARAM,sim_params.rng);
		}
	}
}



void copy_arrays_mcbf(Simulation_Parameters& sim_params, double* j_to,  double* k_to, double* b_to, double *j_from, double *k_from, double* b_from, int destination_index, int source_index){


	std::memcpy(&j_to[2*destination_index], &j_from[2*source_index], 2*NUMBER_OF_SITES*sim_params.total_steps*sizeof(double));
	std::memcpy(&k_to[2*destination_index], &k_from[2*source_index], 2*NUMBER_OF_SITES*sim_params.total_steps*sizeof(double));
	std::memcpy(&b_to[destination_index], &b_from[source_index], NUMBER_OF_SITES*sim_params.total_steps*sizeof(double));
}



double get_change_mcbf(Simulation_Parameters& sim_params){
	int sign;
	double max_change, abs_change, temp_scalar;

  sign       = pow(-1,(int)floor(get_random_double(0,10,sim_params.rng)));
	max_change = MAX_CHANGE_MCBF_INIT * (TAU_INIT/sim_params.tau);
	abs_change = get_random_double(0,max_change, sim_params.rng);

	return sign*abs_change;
}





void change_array_mcbf(double *j_array, double *k_array, double *b_array, int random_row, int random_col, double change, int i, int total_steps){
	int mod = VARIATIONS;


	if(i%mod==ROW-1) change_row_mcbf(j_array,k_array,b_array,random_row,change,true, true,true, 1, 0);
	else if(i%mod==COL-1) change_col_mcbf(total_steps,j_array,k_array,b_array,random_col,change, true, true, true, 1, 0);
	else if(i%mod==ALTROW-1) change_row_mcbf(j_array,k_array,b_array,random_row,change, true, true, true ,2, 0);
	else if(i%mod==ALTCOL-1) change_col_mcbf(total_steps,j_array,k_array,b_array,random_col,change, true, true, true, 2, 0);
	else if(i%mod==ALTROW2-1) change_row_mcbf(j_array,k_array,b_array,random_row,change, true, true, true, 2, 1);
	else if(i%mod==ALTCOL2-1) change_col_mcbf(total_steps,j_array,k_array,b_array,random_col,change, true, true, true, 2, 1);
	else if(i%mod==ROWJ-1) change_row_mcbf(j_array,k_array,b_array,random_row,change, true, false, false, 1, 0);
	else if(i%mod==ROWK-1) change_row_mcbf(j_array,k_array,b_array,random_row,change, false, true, false, 1, 0);
	else if(i%mod==ROWB-1) change_row_mcbf(j_array,k_array,b_array,random_row,change, false, false, true, 1, 0);
	else if(i%mod==COLJ-1) change_col_mcbf(total_steps,j_array,k_array,b_array,random_col,change, true, false, false, 2, 1);
	else if(i%mod==COLK-1) change_col_mcbf(total_steps,j_array,k_array,b_array,random_col,change, false, true, false, 2, 1);
	else if(i%mod==COLB-1) change_col_mcbf(total_steps,j_array,k_array,b_array,random_col,change, false, false, true, 2, 1);
	else if(i%mod==SINGLE-1) change_single_mcbf(j_array,k_array,b_array,random_row,random_col,change);
	else printf("ERROR IN CHANGE_ARRAY_MCBF\n"), exit(0);
}



void change_row_mcbf(double *j_array,double *k_array,double *b_array, int row, double change, bool j, bool k, bool b, int jump, int offset){
	int i;


	if(k) for (i=offset; i<NUMBER_OF_SITES*2; i+=jump){if ((MIN_PARAM < k_array[NUMBER_OF_SITES*2*row+i] + change) &&  (k_array[NUMBER_OF_SITES*2*row+i] + change < MAX_PARAM)) k_array[NUMBER_OF_SITES*2*row+i] += change;}
	if(j) for (i=offset; i<NUMBER_OF_SITES*2; i+=jump){if ((MIN_PARAM < j_array[NUMBER_OF_SITES*2*row+i] + change) && (j_array[NUMBER_OF_SITES*2*row+i] + change < MAX_PARAM)) j_array[NUMBER_OF_SITES*2*row+i] += change;}
	if(b) for (i=offset; i<NUMBER_OF_SITES; i+=jump){if ((MIN_PARAM < b_array[NUMBER_OF_SITES*row+i] + change) && (b_array[NUMBER_OF_SITES*row+i] + change < MAX_PARAM)) b_array[NUMBER_OF_SITES*row+i] += change;}
}



void change_col_mcbf(int total_steps,double *j_array,double *k_array,double *b_array, int col, double change,bool j, bool k, bool b, int jump, int offset){
	int i;


	if(j) for (i=offset; i<total_steps; i+=jump){if ((MIN_PARAM < j_array[NUMBER_OF_SITES*2*i+col] + change) && (j_array[NUMBER_OF_SITES*2*i+col] + change < MAX_PARAM)) j_array[NUMBER_OF_SITES*2*i+col] += change;}
	if(k) for (i=offset; i<total_steps; i+=jump){if ((MIN_PARAM < k_array[NUMBER_OF_SITES*2*i+col] + change) && (k_array[NUMBER_OF_SITES*2*i+col] + change  < MAX_PARAM)) k_array[NUMBER_OF_SITES*2*i+col] += change;}
	if(b) for (i=offset; i<total_steps; i+=jump){if ((MIN_PARAM < b_array[NUMBER_OF_SITES*i+(int)floor(col/2.0)] + change) && (b_array[NUMBER_OF_SITES*i+(int)floor(col/2.0)] + change < MAX_PARAM)) b_array[NUMBER_OF_SITES*i+(int)floor(col/2.0)] += change;}
}



void change_single_mcbf(double *j_array,double *k_array,double *b_array, int row,int col, double change){


	if ((MIN_PARAM < j_array[NUMBER_OF_SITES*2*row+col] + change) && (j_array[NUMBER_OF_SITES*2*row+col] + change < MAX_PARAM)) j_array[NUMBER_OF_SITES*2*row+col] += change;
 	if ((MIN_PARAM < k_array[NUMBER_OF_SITES*2*row+col] + change) && (k_array[NUMBER_OF_SITES*2*row+col] + change < MAX_PARAM)) k_array[NUMBER_OF_SITES*2*row+col] += change;
	if ((MIN_PARAM < b_array[NUMBER_OF_SITES*row+(int)floor(col/2.0)] + change) && (b_array[NUMBER_OF_SITES*row+(int)floor(col/2.0)] + change  < MAX_PARAM)) b_array[NUMBER_OF_SITES*row+(int)floor(col/2.0)] += change;
}




void scale_arrays_mcbf(double* j_array, double* k_array,double* b_array, int total_steps){
	// int i, j;
	//
	//
	// for(i=total_steps;i>0;i--){
	// 	for(j=0;j<NUMBER_OF_SITES*2;j++){
	// 		j_array[NUMBER_OF_SITES*2*(ARRAY_SCALAR*i-1)+j] = j_array[NUMBER_OF_SITES*2*(i-1)+j];
	// 		j_array[NUMBER_OF_SITES*2*(ARRAY_SCALAR*i-2)+j] = j_array[NUMBER_OF_SITES*2*(i-1)+j];
	// 		k_array[NUMBER_OF_SITES*2*(ARRAY_SCALAR*i-1)+j] = k_array[NUMBER_OF_SITES*2*(i-1)+j];
	// 		k_array[NUMBER_OF_SITES*2*(ARRAY_SCALAR*i-2)+j] = k_array[NUMBER_OF_SITES*2*(i-1)+j];
	// 	}
	// 	for(j=0;j<NUMBER_OF_SITES;j++){
	// 		b_array[NUMBER_OF_SITES*(ARRAY_SCALAR*i-1)+j] = 	b_array[NUMBER_OF_SITES*(i-1)+j];
	// 		b_array[NUMBER_OF_SITES*(ARRAY_SCALAR*i-2)+j] = 	b_array[NUMBER_OF_SITES*(i-1)+j];
	// 	}
	// }
}








void binary_search_mcbf(Simulation_Parameters& sim_params){
	// int i;
	// double tau_max = sim_params.tau, tau_min;
	//
	// tau_min = sim_params.tau/TAU_SCALAR_TINY;
	// sim_params.tau = (tau_max + tau_min) / 2.0;
	// sim_params.new_distance = sim_params.old_distance;
	//
	// if(PRINT) printf("\nUsing binary search method to look for optimal ground state...");
	//
	// while((tau_max - tau_min) >BINARY_SEARCH_TAU_LIMIT_MCBF){
	//
	// 	gsl_rng_set(sim_params.rng, 1);
	// 	calc_initial_temp_mcbf(sim_params);
	//
	// 	for(sim_params.seed=1; sim_params.seed<NUM_SEEDS+1; sim_params.seed++){
	// 		gsl_rng_set(sim_params.rng, sim_params.seed);
	// 		mcbf_simulation(sim_params);
	//
	// 		sim_params.E_array_fixed_tau[sim_params.seed-1] = sim_params.best_E;
	// 		copy_arrays_mcbf(sim_params,sim_params.j_best_fixed_tau,sim_params.k_best_fixed_tau,sim_params.b_best_fixed_tau,sim_params.j_best,  sim_params.k_best,  sim_params.b_best,  (sim_params.seed-1)*NUMBER_OF_SITES*sim_params.total_steps, 0);
	// 	}
	// 	for(sim_params.seed=1; sim_params.seed<NUM_SEEDS+1; sim_params.seed++) if(sim_params.E_array_fixed_tau[sim_params.seed-1] < sim_params.best_E) sim_params.best_E = sim_params.E_array_fixed_tau[sim_params.seed-1];
	//
	// 	sim_params.old_distance = sim_params.new_distance;
	// 	sim_params.new_distance = calc_distance(sim_params.intial_state, sim_params.state);
	//
	// 	if(PRINT) print_mc_results(sim_params);
	//
	// 	if(sim_params.new_distance < DISTANCE_LIMIT){
	// 		if(PRINT) printf("\nIn binary_search_mcbf....\nStepping backward....\n");
	// 		tau_max = sim_params.tau;
	// 		sim_params.tau = (tau_max+tau_min)/2.0;
	// 		sim_params.time_step = sim_params.tau/((double) sim_params.total_steps);
	//
	// 	}else{
	// 		if(PRINT) printf("\nIn binary_search_mcbf....\nStepping forward....\n");
	// 		tau_min = sim_params.tau;
	// 		sim_params.tau = (tau_min+tau_max)/2;
	// 		calc_total_steps_mcbf(sim_params);
	// 		sim_params.time_step = sim_params.tau/((double) sim_params.total_steps);
	// 		if(SAVE_DATA) save_mcbf_data_fixed_tau(sim_params);
	// 	}
	// }
	// if(PRINT) printf("\nUpper and lower tau are close enough, exiting binary_search_mcbf\n");
}
