#include <gsl/gsl_rng.h>
#include <ctime>
#include <cstring>

#include "parameters.h"
#include "operations.h"
#include "check.h"
#include "hamiltonian.h"
#include "linear_algebra.h"




void Simulation_Parameters::initialize_lattice(int number_of_occupants, Simulation_Parameters& sim_params){
	const gsl_rng_type * TT;

	num_occupants = number_of_occupants;
	N	      = choose(num_occupants);
	b 	      = new unsigned long long int[N]();
	ham_target    = new double[2*N*N]();
	ham_initial   = new double[2*N*N]();
	init_state   = new double[2*N]();
	target_state  = new double[2*N]();
	jkb_initial   = new double[3]();
	jkb_target    = new double[3]();
	table 	      = new int[num_occupants*N]();
	bonds	        = new int[NUMBER_OF_SITES*NUMBER_OF_SITES]();
	evolved_state_fixed_tau = new double[2*N*NUM_SEEDS]();
	best_evolved_state  = new double[2*N]();

	construct_lattice(lattice);
	assign_bonds(bonds, lattice);
	combinations(N, num_occupants,b,table);
	gsl_rng_env_setup();
	TT = gsl_rng_default;
	rng  = gsl_rng_alloc (TT);
}

void Simulation_Parameters::initialize_hamiltonians(double ji, double ki, double jt, double kt, Simulation_Parameters& sim_params){
	int INCX = 1,INCY = 1;


	j_initial = ji;
	k_initial = ki;
	b_initial = 0;

	j_target = jt;
	k_target = kt;
	b_target = 0;

	jkb_initial[0] = j_initial;
	jkb_initial[1] = k_initial;
	jkb_initial[2] = b_initial;

	jkb_target[0] = j_target;
	jkb_target[1] = k_target;
	jkb_target[2] = b_target;

	construct_device_hamiltonian_uniform(sim_params, ham_initial, jkb_initial);
	construct_device_hamiltonian_uniform(sim_params, ham_target, jkb_target);

	get_ground_state(N, ham_initial, init_state);
	get_ground_state(N, ham_target, target_state);
	ground_E = cost(target_state, N, target_state, ham_target, true);
	initial_E = cost(target_state, N, init_state, ham_target, true);
	init_target_dot_squared = pow(zdotc_(&N, target_state, &INCX, init_state, &INCY),2);

	evolved_target_dot_squared = 0;

	if(CHECK) check_norm(init_state, N);
}



void Simulation_Parameters::init_mcbb_params(Simulation_Parameters& sim_params){
	best_mc_result_fixed_tau      = new double[NUM_SEEDS]();
	j_best_fixed_tau       = new double[NUM_SEEDS*2*NUMBER_OF_BANGS]();
	k_best_fixed_tau       = new double[NUM_SEEDS*2*NUMBER_OF_BANGS]();
	b_best_fixed_tau       = new double[NUM_SEEDS*2*NUMBER_OF_BANGS]();

	j_best		          = new double[2*NUMBER_OF_BANGS]();
	k_best  	          = new double[2*NUMBER_OF_BANGS]();
	b_best 		          = new double[2*NUMBER_OF_BANGS]();


	tau		  = TAU_INIT;
	old_distance      = 1;
	new_distance      = 1;
	sweeps_multiplier = 1;

	start = std::clock();

	if(DIAG){
		int i,j;
		double *filler_ham, *real_ham;
		filler_ham                = new double[2*N*N]();
		real_ham                  = new double[N*N]();
		P_11 		          = new double[2*N*N]();
		P_11_inv                  = new double[2*N*N]();
		D_11                      = new double[N]();
		P_10 		          = new double[2*N*N]();
		P_10_inv                  = new double[2*N*N]();
		D_10                      = new double[N]();
		P_01 		          = new double[2*N*N]();
		P_01_inv                  = new double[2*N*N]();
		D_01                      = new double[N]();

		double jkb_11[3] = {1,1,0};
		double jkb_10[3] = {1,0,0};
		double jkb_01[3] = {0,1,0};
		construct_device_hamiltonian_uniform(sim_params, filler_ham, jkb_11);
		for(i=0;i<N*N;i++) real_ham[i] = filler_ham[i*2];
		diag_hermitian_real_double(N, real_ham,P_11,D_11);
		construct_device_hamiltonian_uniform(sim_params, filler_ham, jkb_10);
		for(i=0;i<N*N;i++) real_ham[i] = filler_ham[i*2];
		diag_hermitian_real_double(N, real_ham,P_10,D_10);
		construct_device_hamiltonian_uniform(sim_params, filler_ham, jkb_01);
		for(i=0;i<N*N;i++) real_ham[i] = filler_ham[i*2];
		diag_hermitian_real_double(N, real_ham,P_01,D_01);
		for(i=0;i<N;i++){
			for(j=0;j<N;j++){
				P_11_inv[i+j*N] = P_11[i*N+j];
				P_10_inv[i+j*N] = P_10[i*N+j];
				P_01_inv[i+j*N] = P_01[i*N+j];
			}
		}
		for(i=0;i<N*N;i++){
			P_11_inv[2*N*N-(i+1)*2] = P_11_inv[N*N-(i+1)];
			P_11_inv[2*N*N-(i+1)*2+1] = 0;
			P_01_inv[2*N*N-(i+1)*2] = P_01_inv[N*N-(i+1)];
			P_01_inv[2*N*N-(i+1)*2+1] = 0;
			P_10_inv[2*N*N-(i+1)*2] = P_10_inv[N*N-(i+1)];
			P_10_inv[2*N*N-(i+1)*2+1] = 0;
			P_11[2*N*N-(i+1)*2] = P_11[N*N-(i+1)];
			P_11[2*N*N-(i+1)*2+1] = 0;
			P_01[2*N*N-(i+1)*2] = P_01[N*N-(i+1)];
			P_01[2*N*N-(i+1)*2+1] = 0;
			P_10[2*N*N-(i+1)*2] = P_10[N*N-(i+1)];
			P_10[2*N*N-(i+1)*2+1] = 0;
		}
		delete[] filler_ham, delete[] real_ham;

	}
}





void Simulation_Parameters::clear_mcbb_params(){
	delete[] best_mc_result_fixed_tau;
	delete[] j_best_fixed_tau;
	delete[] k_best_fixed_tau;
	delete[] b_best_fixed_tau;
	delete[] j_best;
	delete[] k_best;
	delete[] b_best;
	if(DIAG){

		delete[] P_11;
		delete[] P_11_inv;
		delete[] D_11;    
		delete[] P_10; 	 
		delete[] P_10_inv;
		delete[] D_10;   
		delete[] P_01; 	 
		delete[] P_01_inv;
		delete[] D_01;    
	}
}






void Simulation_Parameters::init_mcdb_params(Simulation_Parameters& sim_params){
	max_steps_mcdb    = MAX_STEPS_MCDB;
    saved_states      = new double[2*N*(MAX_STEPS_MCDB+1)]();
	best_mc_result_fixed_tau = new double[NUM_SEEDS]();
	j_best_fixed_tau  = new double[NUM_SEEDS*max_steps_mcdb]();
	k_best_fixed_tau  = new double[NUM_SEEDS*max_steps_mcdb]();
	b_best_fixed_tau  = new double[NUM_SEEDS*max_steps_mcdb]();

	j_best  	  = new double[max_steps_mcdb]();
	k_best  	  = new double[max_steps_mcdb]();
	b_best 		  = new double[max_steps_mcdb]();

	j_best_scaled     = new double[max_steps_mcdb]();
	k_best_scaled 	  = new double[max_steps_mcdb]();
	b_best_scaled	  = new double[max_steps_mcdb]();

	tau                = TAU_INIT;
	total_steps       = MIN_STEPS_MCDB;
	time_step         = tau/total_steps;
	old_distance      = 1;
	new_distance      = 1;
	sweeps_multiplier = 1;

	e10  = new double[2*N*N*(TOTAL_STEP_CHANGES+1)]();
	e01  = new double[2*N*N*(TOTAL_STEP_CHANGES+1)]();
	e11  = new double[2*N*N*(TOTAL_STEP_CHANGES+1)]();

	start = std::clock();
	for(int k=0; k<2*N; k++) saved_states[k] = init_state[k];
	if(MCBB_SECONDARY){
		j_best_secondary      = new double[2*NUMBER_OF_BANGS]();
		k_best_secondary      = new double[2*NUMBER_OF_BANGS]();
		b_best_secondary      = new double[2*NUMBER_OF_BANGS]();
		if(DIAG){
			int i,j;
			double *filler_ham, *real_ham;
			filler_ham                = new double[2*N*N]();
			real_ham                  = new double[N*N]();
			P_11 		          = new double[2*N*N]();
			P_11_inv                  = new double[2*N*N]();
			D_11                      = new double[N]();
			P_10 		          = new double[2*N*N]();
			P_10_inv                  = new double[2*N*N]();
			D_10                      = new double[N]();
			P_01 		          = new double[2*N*N]();
			P_01_inv                  = new double[2*N*N]();
			D_01                      = new double[N]();

			double jkb_11[3] = {1,1,0};
			double jkb_10[3] = {1,0,0};
			double jkb_01[3] = {0,1,0};
			construct_device_hamiltonian_uniform(sim_params, filler_ham, jkb_11);
			for(i=0;i<N*N;i++) real_ham[i] = filler_ham[i*2];
			diag_hermitian_real_double(N, real_ham,P_11,D_11);
			construct_device_hamiltonian_uniform(sim_params, filler_ham, jkb_10);
			for(i=0;i<N*N;i++) real_ham[i] = filler_ham[i*2];
			diag_hermitian_real_double(N, real_ham,P_10,D_10);
			construct_device_hamiltonian_uniform(sim_params, filler_ham, jkb_01);
			for(i=0;i<N*N;i++) real_ham[i] = filler_ham[i*2];
			diag_hermitian_real_double(N, real_ham,P_01,D_01);
			for(i=0;i<N;i++){
				for(j=0;j<N;j++){
					P_11_inv[i+j*N] = P_11[i*N+j];
					P_10_inv[i+j*N] = P_10[i*N+j];
					P_01_inv[i+j*N] = P_01[i*N+j];
				}
			}
			for(i=0;i<N*N;i++){
				P_11_inv[2*N*N-(i+1)*2] = P_11_inv[N*N-(i+1)];
				P_11_inv[2*N*N-(i+1)*2+1] = 0;
				P_01_inv[2*N*N-(i+1)*2] = P_01_inv[N*N-(i+1)];
				P_01_inv[2*N*N-(i+1)*2+1] = 0;
				P_10_inv[2*N*N-(i+1)*2] = P_10_inv[N*N-(i+1)];
				P_10_inv[2*N*N-(i+1)*2+1] = 0;
				P_11[2*N*N-(i+1)*2] = P_11[N*N-(i+1)];
				P_11[2*N*N-(i+1)*2+1] = 0;
				P_01[2*N*N-(i+1)*2] = P_01[N*N-(i+1)];
				P_01[2*N*N-(i+1)*2+1] = 0;
				P_10[2*N*N-(i+1)*2] = P_10[N*N-(i+1)];
				P_10[2*N*N-(i+1)*2+1] = 0;
			}

			delete[] filler_ham, delete[] real_ham;
		}
	}
}




void Simulation_Parameters::clear_mcdb_params(){
	delete[] best_mc_result_fixed_tau;
	delete[] j_best_fixed_tau;
	delete[] k_best_fixed_tau;
	delete[] b_best_fixed_tau;
	delete[] j_best;
	delete[] k_best;
	delete[] b_best;
	delete[] j_best_scaled;
	delete[] k_best_scaled;
	delete[] b_best_scaled;
	delete[] e10;
	delete[] e01;
	delete[] e11;
	if(DIAG){
		delete[] P_11;
		delete[] P_11_inv;
		delete[] D_11;    
		delete[] P_10; 	 
		delete[] P_10_inv;
		delete[] D_10;   
		delete[] P_01; 	 
		delete[] P_01_inv;
		delete[] D_01;    
	}
}




void Simulation_Parameters::init_mcbf_params(){
	best_mc_result_fixed_tau = new double[NUM_SEEDS]();
	j_best_fixed_tau  = new double[NUM_SEEDS*2*NUMBER_OF_SITES*MAX_EVOLVE_STEPS_MCBF]();
	k_best_fixed_tau  = new double[NUM_SEEDS*2*NUMBER_OF_SITES*MAX_EVOLVE_STEPS_MCBF]();
	b_best_fixed_tau  = new double[NUM_SEEDS*NUMBER_OF_SITES*MAX_EVOLVE_STEPS_MCBF]();

	j_best            = new double[2*NUMBER_OF_SITES*MAX_EVOLVE_STEPS_MCBF]();
	k_best            = new double[2*NUMBER_OF_SITES*MAX_EVOLVE_STEPS_MCBF]();
	b_best            = new double[NUMBER_OF_SITES*MAX_EVOLVE_STEPS_MCBF]();

	tau               = TAU_INIT;
	total_steps       = TOTAL_STEPS_INIT_MCBF;
	time_step         = tau/total_steps;
	old_distance      = 1;
	new_distance      = 1;
	sweeps_multiplier = 1;

	start = std::clock();


}




void Simulation_Parameters::clear_mcbf_params(){
	delete[] best_mc_result_fixed_tau;
	delete[] j_best_fixed_tau;
	delete[] k_best_fixed_tau;
	delete[] b_best_fixed_tau;
	delete[] j_best;
	delete[] k_best;
	delete[] b_best;
}
