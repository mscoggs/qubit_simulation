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


#include "adiabatic.h"
#include "check.h"
#include "export_data.h"
#include "hamiltonian.h"
#include "mcbb.h"
#include "mcbf.h"
#include "operations.h"
#include "parameters.h"
#include "print.h"
#include "linear_algebra.h"


/* 
g++ -c main.cpp adiabatic.cpp check.cpp export_data.cpp hamiltonian.cpp linear_algebra.cpp mcbb.cpp mcbf.cpp operations.cpp print.cpp parameters.cpp
g++ main.o adiabatic.o check.o export_data.o hamiltonian.o linear_algebra.o mcbb.o mcbf.o operations.o print.o parameters.o -lgsl -llapack -lblas
./a.out g_initial(double) f_initial(double) g_target(double) f_target(double) 

OR

g++ -o main main.cpp adiabatic.cpp check.cpp export_data.cpp hamiltonian.cpp linear_algebra.cpp mcbb.cpp mcbf.cpp operations.cpp print.cpp parameters.cpp -lgsl -llapack -lblas -std=gnu++11
./main g_initial(double) f_initial(double) g_target(double) f_target(double) 
*/


main (int argc, char *argv[]){
	int gt_i, ft_i, gi_i, fi_i, occupant_index;
	Simulation_Parameters sim_params;
	double g_init[1] = {atof(argv[1])};//{0.1,1};
	double f_init[1] = {atof(argv[2])};//{0.1,1};
	double g_targ[1] = {atof(argv[3])};//{0.5,0.1,1};
	double f_targ[1] = {atof(argv[4])};//{0.5,0.1,1};
	int num_occupants_array[1]= {2};
	
	

	for(occupant_index=0;occupant_index<sizeof(num_occupants_array)/sizeof(int);occupant_index++){
		sim_params.initialize_parameters(num_occupants_array[occupant_index]);

		for(gi_i=0;gi_i<(sizeof(g_init)/sizeof(double));gi_i++){
			for(fi_i=0;fi_i<(sizeof(f_init)/sizeof(double));fi_i++){
				for(gt_i=0;gt_i<(sizeof(g_targ)/sizeof(double));gt_i++){
					for(ft_i=0;ft_i<(sizeof(f_targ)/sizeof(double));ft_i++){
						sim_params.g_initial = g_init[gi_i];
						sim_params.f_initial = f_init[fi_i];
						sim_params.g_target = g_targ[gt_i];
						sim_params.f_target = f_targ[ft_i];
						//if(sim_params.g_initial == sim_params.g_target && sim_params.f_initial == sim_params.f_target) continue;
						sim_params.initialize_target_and_initial_hamiltonian(sim_params);

						if(MCBF) mcbf_method(sim_params);
						if(MCBB) mcbb_method(sim_params);
						if(ADIABATIC)  adiabatic_method(sim_params);
					}}}}
		delete[] sim_params.b, delete[] sim_params.ham_initial, delete[] sim_params.ham_target, delete[] sim_params.table, delete[] sim_params.start_state;
	}
	exit(0);
}
