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

#include "check.h"
#include "write_data.h"
#include "hamiltonian.h"
#include "operations.h"
#include "parameters.h"
#include "print.h"
#include "linear_algebra.h"



main (int argc, char *argv[]){
	double k, j, *evals, *v_diag, *ham_real, gap, *k_array, *j_array, *gap_array;
	int L, num_occupants, i, index =0;
	L = atoi(argv[1]);
	num_occupants = atoi(argv[2]);
	Simulation_Parameters sim_params;
	sim_params.initialize_simluation(num_occupants, 0,0, 0, 0, sim_params);

	k_array = new double[L*L]();
	j_array = new double[L*L]();
	gap_array = new double[L*L]();
	evals   = new double[sim_params.N]();
	printf("N: %i\n\n\n", sim_params.N);
	v_diag   = new double[sim_params.N*sim_params.N]();
	ham_real   = new double[sim_params.N*sim_params.N]();

	for(k=0.0; k<1; k+=1.0/L){
		for(j=0.0; j<1; j+=1.0/L){
			printf("j: %f, k: %f\n", j, k);

			sim_params.jkb_initial[0] = j;
			sim_params.jkb_initial[1] = (1-j)*k;
			sim_params.jkb_initial[2] = (1-j)*(1-k);

			construct_device_hamiltonian_uniform(sim_params, sim_params.ham_initial, sim_params.jkb_initial);
			for(i=0; i<sim_params.N*sim_params.N; i++) ham_real[i] = sim_params.ham_initial[i*2];

			diag_hermitian_real_double(sim_params.N, ham_real, v_diag, evals);

			gap = evals[1] - evals[0];
			k_array[index] = k;
			j_array[index] = j;
			gap_array[index] = gap;
			index ++;
		}
	}
	save_gap_data(sim_params, j_array, k_array, gap_array, L);
}
