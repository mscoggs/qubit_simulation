#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <string>
#include <algorithm>
#include <ctime>
#include <cstring>

#include "adiabatic.h"
#include "check.h"
#include "export_data.h"
#include "hamiltonian.h"
#include "operations.h"
#include "parameters.h"
#include "print.h"
#include "linear_algebra.h"
#include "operations.h"
/*
g++ -o main main.cpp adiabatic.cpp check.cpp export_data.cpp hamiltonian.cpp linear_algebra.cpp mcbb.cpp mcbf.cpp operations.cpp print.cpp parameters.cpp -lgsl -llapack -lblas -std=gnu++11
./main j_initial(double) k_initial(double) j_target(double) k_target(double)

OR

make compile
make run
*/
double find_next_time(double time, double tau,double* j_array, double* k_array,double*  b_array,int*  jkb_index, double*jkb);


main (int argc, char *argv[]){
	Simulation_Parameters sim_params;
	int num_occupants, *jkb_index, i, j;
	double ji,ki,jt,kt;
	double *hamiltonian,*ham_t_i, *ham_real,*exp_matrix,*e_vals, *v_diag, *jkb_vals, time=0, time2=0, time_step;

	num_occupants = 3;
	ji = 0.05;
	ki = 0.05;
	jt = 0.5;
	kt = 0.95;
	sim_params.initialize_lattice(num_occupants,sim_params);
	sim_params.initialize_hamiltonians(ji,ki,jt,kt, sim_params);
	sim_params.init_mcbb_params();
	sim_params.state = new double[2*sim_params.N]();
	sim_params.tau = 0.5259985;


	std::memcpy(sim_params.state,sim_params.start_state, 2*sim_params.N*sizeof(double));
	double j_array[6] = {0,0,0.053341200000000005, 0.207361, 0.2827445, 0.42718500000000004};
	double k_array[6] = {0,0,0,0,0, 0.51};
	double b_array[6] = {0,0,0,0,0, 0};

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
	sim_params.best_E = cost(sim_params.N, sim_params.state, sim_params.ham_target);
	printf("BEST_E: %f\n\n", sim_params.best_E);
	delete[] hamiltonian, delete[] ham_t_i, delete[] exp_matrix, delete[] e_vals, delete[] v_diag, delete[] ham_real, delete[] jkb_vals, delete[] jkb_index;
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
