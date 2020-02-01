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
#include "export_data.h"
#include "hamiltonian.h"
#include "operations.h"
#include "parameters.h"
#include "print.h"
#include "linear_algebra.h"

/*
g++ -o main main.cpp adiabatic.cpp check.cpp export_data.cpp hamiltonian.cpp linear_algebra.cpp mcbb.cpp mcbf.cpp operations.cpp print.cpp parameters.cpp -lgsl -llapack -lblas -std=gnu++11
./main g_initial(double) f_initial(double) g_target(double) f_target(double) 
*/

main (int argc, char *argv[]){
	double f, g, *evals, *v_diag, *ham_real, gap, *f_array, *g_array, *gap_array;
	int L, num_occupants, i, index =0;
	L = atoi(argv[1]);
	num_occupants = atoi(argv[2]);
	Simulation_Parameters sim_params;
	sim_params.initialize_parameters(num_occupants);

	f_array = new double[L*L]();
	g_array = new double[L*L]();
	gap_array = new double[L*L]();
	evals   = new double[sim_params.N]();
	printf("N: %i\n\n\n", sim_params.N);
	v_diag   = new double[sim_params.N*sim_params.N]();
	ham_real   = new double[sim_params.N*sim_params.N]();

	for(f=0.0; f<1; f+=1.0/L){
		for(g=0.0; g<1; g+=1.0/L){
			printf("g: %f, f: %f\n", g, f);
			sim_params.jkb_initial[0] = g;
			sim_params.jkb_initial[1] = (1-g)*f;
			sim_params.jkb_initial[2] = (1-g)*(1-f);

			construct_device_hamiltonian_uniform(sim_params, sim_params.ham_initial, sim_params.jkb_initial);
			for(i=0; i<sim_params.N*sim_params.N; i++) ham_real[i] = sim_params.ham_initial[i*2];

			diag_hermitian_real_double(sim_params.N, ham_real, v_diag, evals);
			
			gap = evals[1] - evals[0];
			f_array[index] = f;
			g_array[index] = g;
			gap_array[index] = gap;
			index ++;
		}
	}	
	save_gap_data(sim_params, f_array, g_array, gap_array, L);
}
