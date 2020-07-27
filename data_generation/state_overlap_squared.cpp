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
#include <iostream>

#include "check.h"
#include "export_data.h"
#include "hamiltonian.h"
#include "operations.h"
#include "parameters.h"
#include "print.h"
#include "linear_algebra.h"
#include "operations.h"



main (int argc, char *argv[]){
	Simulation_Parameters sim_params;
	int a,b,c,d,e;
	std::ofstream file;

	file.open("../data/overlap.txt");
	file << "num_occupants[a]  ji[b] ki[c]  jt[d]  kt[e]  sim_params.init_target_dot_squared\n";

	// int num_occupants[4] = {2,3,4,7};
	// double ji[2] = {0.050, 0.95};
	// double ki[2] = {0.050, 0.95};
	// double jt[13] = {0.050, 0.125, 0.20, 0.275, 0.35, 0.425, 0.50, 0.575, 0.65, 0.725, 0.8, 0.875, 0.95};
	// double kt[13] = {0.050, 0.125, 0.20, 0.275, 0.35, 0.425, 0.50, 0.575, 0.65, 0.725, 0.8, 0.875, 0.95};
	// for(a=0;a<4;a++){
	// 	for(b=0;b<2;b++){
	// 		for(c=0;c<2;c++){
	// 			for(d=0;d<13;d++){
	// 				for(e=0;e<13;e++){

	int num_occupants[2] = {4,7};
	double ji[2] = {0.050, 0.95};
	double ki[2] = {0.050, 0.95};
	double jt[13] = {0.050, 0.125, 0.20, 0.275, 0.35, 0.425, 0.50, 0.575, 0.65, 0.725, 0.8, 0.875, 0.95};
	double kt[13] = {0.050, 0.125, 0.20, 0.275, 0.35, 0.425, 0.50, 0.575, 0.65, 0.725, 0.8, 0.875, 0.95};
	for(a=0;a<2;a++){
		for(b=0;b<2;b++){
			for(c=0;c<2;c++){
				for(d=0;d<13;d++){
					for(e=0;e<13;e++){
						sim_params.initialize_lattice(num_occupants[a],sim_params);
						sim_params.initialize_hamiltonians(ji[b],ki[c],jt[d],kt[e], sim_params);
						std::string jis = std::to_string(ji[b]);
						jis.erase(jis.find_last_not_of('0') + 2, std::string::npos);
						std::string kis = std::to_string(ki[c]);
						kis.erase(kis.find_last_not_of('0') + 2, std::string::npos);
						std::string jts = std::to_string(jt[d]);
						jts.erase(jts.find_last_not_of('0') + 2, std::string::npos);
						std::string kts = std::to_string(kt[e]);
						kts.erase(kts.find_last_not_of('0') + 2, std::string::npos);

						file << num_occupants[a] << " " << jis << " " << kis << " " << jts << " " << kts << " " << sim_params.init_target_dot_squared <<"\n";
					}
				}
			}
		}
	}
}
