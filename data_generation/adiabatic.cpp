#include <gsl/gsl_rng.h>
#include <string.h>

#include "adiabatic.h"
#include "check.h"
#include "export_data.h"
#include "hamiltonian.h"
#include "linear_algebra.h"
#include "operations.h"
#include "parameters.h"
#include "print.h"


void adiabatic_method(Simulation_Parameters& sim_params){
	sim_params.tau_array    = new double[MAX_TAU_STEPS_ADIA]();
	sim_params.best_E_array = new double[MAX_TAU_STEPS_ADIA]();
	sim_params.state        = new double[2*sim_params.N]();

	memcpy(sim_params.state,  sim_params.start_state,  2*sim_params.N*sizeof(double));

	sim_params.tau         = TAU_INIT_ADIA;
	sim_params.time_step   = TIME_STEP_ADIA;
	sim_params.total_steps = floor(int(TAU_INIT_ADIA/double(TIME_STEP_ADIA)));
	sim_params.index=0;

	if(PRINT) print_adiabatic_info(sim_params);

	while(sim_params.tau<MAX_TAU_ADIA){
		memcpy(sim_params.state,sim_params.start_state, 2*sim_params.N*sizeof(double));//resetting state
		evolve_adiabatic(sim_params);
		sim_params.best_E_array[sim_params.index] = cost(sim_params.N, sim_params.state, sim_params.ham_target);
		sim_params.tau_array[sim_params.index] = sim_params.tau;

		printf("       Tau: %5.2f  ||  Expectation-Value:  %7.4f  ||  Target: %f \n", sim_params.tau, sim_params.best_E_array[sim_params.index],sim_params.ground_E);

		if((sim_params.best_E_array[sim_params.index] - sim_params.ground_E) < DIFFERENCE_LIMIT_ADIA){
			printf("GROUND STATE REACHED. STOPPING THE ADIABATIC EVOLUTION\n\n");
		  break;
		}

		sim_params.index++;
		sim_params.tau         = sim_params.tau+TAU_INIT_ADIA;//I*TAU_SCALAR_ADIA;
		sim_params.total_steps = floor(sim_params.tau/sim_params.time_step);
	}
	if (ADIABATIC_DATA) save_adiabatic_data(sim_params);
	delete[] sim_params.tau_array, delete[] sim_params.best_E_array, delete[] sim_params.state;
}



void evolve_adiabatic(Simulation_Parameters& sim_params){
	int i,j;
	double *hamiltonian,*ham_t_i, *ham_real,*exp_matrix,*v_diag, *e_vals,*jkb,g,f;

	hamiltonian        = new double[2*sim_params.N*sim_params.N]();
	ham_t_i    = new double[2*sim_params.N*sim_params.N]();
	ham_real   = new double[sim_params.N*sim_params.N]();
	exp_matrix = new double[2*sim_params.N*sim_params.N]();
	v_diag     = new double[sim_params.N*sim_params.N]();
	e_vals     = new double[sim_params.N]();
	jkb        = new double[3]();

	for (i=1; i<sim_params.total_steps+1;i++){

		g = sim_params.g_initial + (sim_params.g_target-sim_params.g_initial)*((i*1.0)/sim_params.total_steps);
		f = sim_params.f_initial + (sim_params.f_target-sim_params.f_initial)*(i*1.0/sim_params.total_steps);
		jkb[0] = g, jkb[1] = (1-g)*f, jkb[2] = (1-f)*(1-g);
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
	delete[] hamiltonian, delete[] ham_t_i, delete[] ham_real, delete[] exp_matrix, delete[] v_diag, delete[] e_vals, delete[] jkb;
}
