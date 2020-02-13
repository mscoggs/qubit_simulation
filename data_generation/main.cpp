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

#include "adiabatic.h"
#include "check.h"
#include "export_data.h"
#include "hamiltonian.h"
#include "mcbb.h"
#include "mcdb.h"
#include "mcbf.h"
#include "operations.h"
#include "parameters.h"
#include "print.h"
#include "linear_algebra.h"

/*
g++ -o main main.cpp adiabatic.cpp check.cpp export_data.cpp hamiltonian.cpp linear_algebra.cpp mcbb.cpp mcbf.cpp operations.cpp print.cpp parameters.cpp -lgsl -llapack -lblas -std=gnu++11
./main g_initial(double) f_initial(double) g_target(double) f_target(double)
*/


main (int argc, char *argv[]){
	Simulation_Parameters sim_params;
	int num_occupants;
	double gi, fi, ft, gt;
	num_occupants = atoi(argv[1]);
	gi = atof(argv[2]);
	fi = atof(argv[3]);
	ft = atof(argv[4]);
	gt = atof(argv[5]);



	sim_params.initialize_simluation(num_occupants, gi, fi, gt, ft, sim_params);
	if(ADIABATIC)  adiabatic_method(sim_params);
	if(MCDB) mcdb_method(sim_params);
	if(MCBB) mcbb_method(sim_params);
	if(MCBF) mcbf_method(sim_params);
}
