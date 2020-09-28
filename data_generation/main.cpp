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
#include "write_data.h"
#include "hamiltonian.h"
#include "mcbb.h"
#include "mcdb.h"
#include "mcbf.h"
#include "operations.h"
#include "parameters.h"
#include "print.h"
#include "linear_algebra.h"



main (int argc, char *argv[]){

	Simulation_Parameters sim_params;
	int num_occupants;
	double ji,ki,jt,kt;
	num_occupants = atoi(argv[1]);
	ji = atof(argv[2]);
	ki = atof(argv[3]);
	jt = atof(argv[4]);
	kt = atof(argv[5]);
	printf("TEST\n");

	sim_params.initialize_lattice(num_occupants,sim_params);
	printf("TEST\n");

	sim_params.initialize_hamiltonians(ji,ki,jt,kt, sim_params);
	printf("TEST\n");

	if(ADIA)  adiabatic_method(sim_params);
	if(MCDB) mcdb_method(sim_params);
	if(MCBB) mcbb_method(sim_params);
	if(MCBF) mcbf_method(sim_params);
}
