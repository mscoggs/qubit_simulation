#include <gsl/gsl_rng.h>
#include <string.h>

#include "adiabatic.h"
#include "check.h"
#include "write_data.h"
#include "hamiltonian.h"
#include "linear_algebra.h"
#include "operations.h"
#include "parameters.h"
#include "print.h"


void adiabatic_method(Simulation_Parameters& sim_params){


	if(check_commutator(sim_params.N, sim_params.ham_initial, sim_params.ham_target) || (sim_params.init_target_dot_squared > INIT_OVERLAP_LIMIT && USE_ENERGY_DISTANCE) || (1-sim_params.init_target_dot_squared < DISTANCE_LIMIT && !USE_ENERGY_DISTANCE)){
		sim_params.tau = 0.0, sim_params.new_distance = 0.0, sim_params.best_mc_result = 0.0;
		if(SAVE_DATA) save_adiabatic_data(sim_params);
		return;
	}

	sim_params.state        = new double[2*sim_params.N]();
	sim_params.start        = std::clock();
	memcpy(sim_params.state,  sim_params.init_state,  2*sim_params.N*sizeof(double));

	sim_params.tau         = TAU_INIT;
	sim_params.total_steps = TOTAL_STEPS_ADIA;//floor(int(TAU_INIT/double(TIME_STEP_ADIA)));
	sim_params.time_step   = sim_params.tau/sim_params.total_steps;//TIME_STEP_ADIA;
	sim_params.new_distance = 1;
	sim_params.old_distance = 1;

	if(PRINT) print_adiabatic_info(sim_params);

	while(sim_params.tau<MAX_TAU_ADIA){
		memcpy(sim_params.state,sim_params.init_state, 2*sim_params.N*sizeof(double));//resetting state
		evolve_adiabatic(sim_params);
		sim_params.best_mc_result = cost(sim_params.target_state, sim_params.N, sim_params.state, sim_params.ham_target);
  	sim_params.old_distance = sim_params.new_distance;
		sim_params.new_distance = calc_distance(sim_params);
		printf("test\n");
		if(PRINT) printf("       Tau: %5.2f  ||  Cost Output:  %7.4f  ||  Distance: %f \n", sim_params.tau, sim_params.best_mc_result,sim_params.new_distance);

		if (SAVE_DATA){
			printf("test2\n");
			sim_params.duration = (std::clock() - sim_params.start)/(double) CLOCKS_PER_SEC;
			save_adiabatic_data(sim_params);
		}
		if((sim_params.new_distance < DISTANCE_LIMIT && USE_ENERGY_DISTANCE) || (1-sim_params.evolved_target_dot_squared < DISTANCE_LIMIT && !USE_ENERGY_DISTANCE)) break;

		sim_params.tau         = sim_params.tau*TAU_SCALAR_ADIA;
		sim_params.time_step   = sim_params.tau/sim_params.total_steps;//TIME_STEP_ADIA;
	}
	delete[] sim_params.state;
}



void evolve_adiabatic(Simulation_Parameters& sim_params){
	int i,j, INCX = 1,INCY = 1;
	double *hamiltonian,*ham_t_i, *ham_real,*exp_matrix,*v_diag, *e_vals,*jkb,j_,k_;

	hamiltonian        = new double[2*sim_params.N*sim_params.N]();
	ham_t_i    = new double[2*sim_params.N*sim_params.N]();
	ham_real   = new double[sim_params.N*sim_params.N]();
	exp_matrix = new double[2*sim_params.N*sim_params.N]();
	v_diag     = new double[sim_params.N*sim_params.N]();
	e_vals     = new double[sim_params.N]();
	jkb        = new double[3]();

	for (i=1; i<sim_params.total_steps+1;i++){

		j_ = sim_params.j_initial + (sim_params.j_target-sim_params.j_initial)*((i*1.0)/sim_params.total_steps);
		k_ = sim_params.k_initial + (sim_params.k_target-sim_params.k_initial)*(i*1.0/sim_params.total_steps);
		jkb[0] = j_, jkb[1] = k_, jkb[2] = 0;
		construct_device_hamiltonian_uniform(sim_params,hamiltonian, jkb);

		if(DIAG){

			for (j=0; j<sim_params.N*sim_params.N; j++) ham_real[j]=0.0, v_diag[j]=0.0, ham_real[j] = hamiltonian[2*j];//converting an all-real-valued complex matrix into just real matrix
			for (j=0; j<sim_params.N; j++) e_vals[j]=0.0;

			diag_hermitian_real_double(sim_params.N, ham_real,v_diag, e_vals);
			exp_diaganolized_real_matrix(hamiltonian, v_diag, e_vals, sim_params.N, sim_params.time_step);//This function exponentiates e_vals to e^(-i*time_step*e_vals)

			if(CHECK) check_unitary(hamiltonian, sim_params.N);
			matrix_vector_mult(hamiltonian,sim_params.state, sim_params.N);

		}else{

			for (j=0; j<sim_params.N*sim_params.N*2; j++) exp_matrix[j] = 0.0, ham_t_i[j]=0.0;
			for (j=0; j<sim_params.N*sim_params.N; j++){ ham_t_i[2*j] = (hamiltonian[2*j+1]*sim_params.time_step), ham_t_i[2*j+1] = (hamiltonian[2*j]*-sim_params.time_step);} //multiplying by -i*dt for the Pade approximation
			exp_complex_double_matrix_pade(sim_params.N, ham_t_i, exp_matrix);

			if(CHECK) check_unitary(exp_matrix, sim_params.N);
			matrix_vector_mult(exp_matrix, sim_params.state, sim_params.N);
		}

		if(CHECK) check_norm(sim_params.state, sim_params.N);

	}

	memcpy(sim_params.best_evolved_state,sim_params.state, 2*sim_params.N*sizeof(double));
	sim_params.evolved_target_dot_squared = complex_dot_squared(sim_params.N*2, sim_params.target_state, sim_params.best_evolved_state);

	delete[] hamiltonian, delete[] ham_t_i, delete[] ham_real, delete[] exp_matrix, delete[] v_diag, delete[] e_vals, delete[] jkb;
}
