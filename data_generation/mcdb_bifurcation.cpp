#include <gsl/gsl_rng.h>
#include <cstring>
#include <algorithm>    // std::min
#include <fstream>
#include <string>
#include <iostream>

#include "write_data.h"
#include "mcdb.h"
#include "mcbb.h"
#include "check.h"
#include "hamiltonian.h"
#include "linear_algebra.h"
#include "operations.h"
#include "parameters.h"
#include "print.h"


void mcdb_method_bifurcation(Simulation_Parameters& sim_params, double time);
void save_mcdb_data_bifurcation(Simulation_Parameters& sim_params);








main (int argc, char *argv[]){

	Simulation_Parameters sim_params;
	int num_occupants,nx;
	double ji,ki,jt,kt;
	num_occupants = atoi(argv[1]);
	ji = atof(argv[2]);
	ki = atof(argv[3]);
	jt = atof(argv[4]);
	kt = atof(argv[5]);
	sim_params.initialize_lattice(num_occupants,sim_params);
	sim_params.initialize_hamiltonians(ji,ki,jt,kt, sim_params);
	double time = atof(argv[6]);
	mcdb_method_bifurcation(sim_params, time);
}


void mcdb_method_bifurcation(Simulation_Parameters& sim_params, double time){
	int i;


	sim_params.init_mcdb_params(sim_params);
	if(check_commutator(sim_params.N, sim_params.ham_initial, sim_params.ham_target) || (1-sim_params.init_target_dot_squared < DISTANCE_LIMIT)){
		sim_params.tau = 0.0, sim_params.new_distance = 0.0, sim_params.best_mc_result = 0.0;
		if(SAVE_DATA) save_mcdb_data_bifurcation(sim_params);
		sim_params.clear_mcdb_params();
		return;
	}

	sim_params.tau = time;
	sim_params.total_steps = MIN_STEPS_MCDB;
	sim_params.time_step = sim_params.tau/sim_params.total_steps;
	pre_exponentiate(sim_params);

	while(sim_params.total_steps <= MAX_STEPS_MCDB){

		gsl_rng_set(sim_params.rng, 0);
		get_sweeps_mcdb(sim_params);
		calc_initial_temp_mcdb(sim_params);

		for(sim_params.seed=1; sim_params.seed<NUM_SEEDS+1; sim_params.seed++){
			gsl_rng_set(sim_params.rng, sim_params.seed);
			mcdb_simulation(sim_params);

			sim_params.best_mc_result_fixed_tau[sim_params.seed-1] = sim_params.best_mc_result;
			std::memcpy(&sim_params.evolved_state_fixed_tau[(sim_params.seed-1)*2*sim_params.N],sim_params.best_evolved_state, 2*sim_params.N*sizeof(double));
			copy_arrays_mcdb(sim_params,sim_params.j_best_fixed_tau,sim_params.k_best_fixed_tau,sim_params.b_best_fixed_tau,sim_params.j_best,  sim_params.k_best,  sim_params.b_best,  (sim_params.seed-1)*sim_params.total_steps, 0);
		}
		get_best_seed(sim_params);

		for(sim_params.seed=1; sim_params.seed<NUM_SEEDS+1; sim_params.seed++) if(sim_params.best_mc_result_fixed_tau[sim_params.seed-1] <= sim_params.best_mc_result) copy_arrays_mcdb(sim_params,sim_params.j_best_scaled,  sim_params.k_best_scaled,  sim_params.b_best_scaled,sim_params.j_best_fixed_tau,sim_params.k_best_fixed_tau,sim_params.b_best_fixed_tau,  0,(sim_params.seed-1)*sim_params.total_steps);
		if(sim_params.total_steps == MAX_STEPS_MCDB) break;
		sim_params.total_steps = sim_params.total_steps*2;
		sim_params.time_step = sim_params.tau/sim_params.total_steps;
		if(sim_params.total_steps <= MAX_STEPS_MCDB) scale_best_arrays_mcdb(sim_params, sim_params.j_best_scaled,sim_params.k_best_scaled,sim_params.b_best_scaled);
	}
	if(MCBB_SECONDARY){
		sim_params.best_mc_result_non_secondary = sim_params.best_mc_result;
		convert_mcdb_to_mcbb(sim_params,sim_params.j_best_secondary,  sim_params.k_best_secondary,  sim_params.b_best_secondary,sim_params.j_best_scaled,  sim_params.k_best_scaled,  sim_params.b_best_scaled);
		mcbb_secondary_simulation(sim_params);
	}
	if(PRINT) print_mc_results(sim_params);
	if(SAVE_DATA) sim_params.duration = (std::clock() - sim_params.start)/(double) CLOCKS_PER_SEC, save_mcdb_data_bifurcation(sim_params);
sim_params.clear_mcdb_params();
}






void save_mcdb_data_bifurcation(Simulation_Parameters& sim_params){
	int i, j;
	std::ofstream file;
	std::string type = "MCDB";
	std::string uni, pbc, ji, ki, jt, kt, ri, rt;

	ji = std::to_string(sim_params.j_initial);
	ji.erase(ji.find_last_not_of('0') + 2, std::string::npos);
	ki = std::to_string(sim_params.k_initial);
	ki.erase(ki.find_last_not_of('0') + 2, std::string::npos);
	jt = std::to_string(sim_params.j_target);
	jt.erase(jt.find_last_not_of('0') + 2, std::string::npos);
	kt = std::to_string(sim_params.k_target);
	kt.erase(kt.find_last_not_of('0') + 2, std::string::npos);

	ri = std::to_string(sim_params.j_initial/sim_params.k_initial);
	rt = std::to_string(sim_params.j_target/sim_params.k_target);

	if (PERIODIC) pbc = 't';
	else pbc = 'f';
	if (UNIFORM_SITES) uni = 't';
	else uni = 'f';

	//std::string dir = "../data/" + std::to_string(NX) + "x" +std::to_string(NY) + "/" + std::to_string(sim_params.num_occupants) ;
	//std::string file_name = "_occupants/" + type + "___PBC="+pbc+"_UNI="+uni+"_DD="+std::to_string(DEVICE_DIMENSION)+"___ji=" + ji + "_ki=" + ki +"_jt=" + jt +  "_kt="+ kt +".txt";
	std::string dir = "../data_bifurcation/" + std::to_string(NX) + "x" +std::to_string(NY) + "/" + std::to_string(sim_params.num_occupants) ;
	std::string file_name = "_occupants/" + type + "___PBC="+pbc+"_UNI="+uni+"_DD="+std::to_string(DEVICE_DIMENSION)+"___ri=" + ri +  "_rt="+ rt +".txt";
	std::string path = dir+file_name;

	file.open(path);
	file << "START_PARAMETERS\n";
	file << "MAX_STEPS_MCDB =               " <<  MAX_STEPS_MCDB << "\n";
	file << "MIN_STEPS_MCDB =               " <<  MIN_STEPS_MCDB << "\n";
	file << "STEPS_CRUNCH_MCDB =            " <<  STEPS_CRUNCH_MCDB << "\n";
	file << "SWEEPS_MCDB =                  " <<  SWEEPS_MCDB << "\n";
	file << "TOTAL_STEP_CHANGES =                  " <<  TOTAL_STEP_CHANGES << "\n";
	file << "MCBB_SECONDARY =                  " <<  MCBB_SECONDARY << "\n";
	file << "SWEEPS_MCBB_SECONDARY =                  " <<  SWEEPS_MCBB_SECONDARY << "\n";
	file << "NUMBER_OF_BANGS =                  " <<  NUMBER_OF_BANGS << "\n";
	file << "GROUND_E =                " <<  sim_params.ground_E << "\n";
	file << "INITIAL_E =               " <<  sim_params.initial_E << "\n";
	file << "j_initial =               " << sim_params.j_initial << "\n";
	file << "k_initial =               " << sim_params.k_initial << "\n";
	file << "b_initial =               " << sim_params.b_initial << "\n";
	file << "j_target =                " << sim_params.j_target << "\n";
	file << "k_target =                " << sim_params.k_target << "\n";
	file << "b_target =                " << sim_params.b_target << "\n";
	file << "init_target_dot_squared = " << sim_params.init_target_dot_squared << "\n";
	file << "init_state = [";
	for(i=0;i<2*sim_params.N;i++) file << sim_params.init_state[i] << ", ";
	file << "]\n";
	file << "target_state = [";
	for(i=0;i<2*sim_params.N;i++) file << sim_params.target_state[i] << ", ";
	file << "]\n";
	file << "USE_ENERGY_DISTANCE =     " << std::boolalpha << USE_ENERGY_DISTANCE << "\n";
	file << "PERIODIC =                " << std::boolalpha << PERIODIC << "\n";
	file << "UNIFORM_SITES =           " << std::boolalpha << UNIFORM_SITES << "\n";
	file << "DEVICE_DIMENSION =        " << DEVICE_DIMENSION << "\n";
	file << "MAX_PARAM =               " << MAX_PARAM << "\n";
	file << "MIN_PARAM =               " << MIN_PARAM << "\n";
	file << "MIN_PARAM =               " << MIN_PARAM << "\n";
	file << "DIAG =                    " <<  std::boolalpha << DIAG << "\n";
	file << "NUM_SEEDS =               " <<  NUM_SEEDS << "\n";
	file << "NUMBER_OF_BANGS =               " <<  NUMBER_OF_BANGS << "\n";
	file << "DISTANCE_LIMIT =          " <<  DISTANCE_LIMIT << "\n";
	file << "TAU_SCALAR =              " <<  TAU_SCALAR << "\n";
	file << "TAU_SCALAR_TINY =         " <<  TAU_SCALAR_TINY << "\n";
	file << "TAU_SCALAR_BIG =          " <<  TAU_SCALAR_BIG << "\n";
	file << "ACCEPTANCE_PROB =         " <<  ACCEPTANCE_PROB << "\n";
	file << "TEMP_EXP_DECAY =          " <<  TEMP_EXP_DECAY << "\n";
	file << "MIN_TEMP_FRACTION =       " <<  MIN_TEMP_FRACTION << "\n";
	file << "TEMP_DECAY_ITERATIONS =   " <<  TEMP_DECAY_ITERATIONS  << "\n";
	file << "ZERO_TEMP_ITREATIONS =    " <<  ZERO_TEMP_ITERATIONS  << "\n";
	file << "RANDOM_STATES =           " <<  RANDOM_STATES  << "\n";
	file << "END_PARAMETERS\n";
	file << "##############################################################################\n";
	file << "clock_duration  =         " << sim_params.duration  << "\n";
	file << "tau =                     " << sim_params.tau   << "\n";
	file << "total_steps =             " << sim_params.total_steps << "\n";
	file << "time_step =               " << sim_params.time_step << "\n";
	file << "best_mc_result =          " << sim_params.best_mc_result << "\n";
	file << "best_mc_result_non_secondary =          " << sim_params.best_mc_result_non_secondary << "\n";
	file << "evolved_target_dot_squared =  " << sim_params.evolved_target_dot_squared << "\n";
	file << "distance =          " << sim_params.new_distance << "\n";

	file << "j_protocol =  [";
	for(i=0;i<NUM_SEEDS;i++){
		file << "[";
		for(j=0;j<sim_params.total_steps;j++) file << sim_params.j_best_fixed_tau[i*sim_params.total_steps + j] << ", ";
		file << "],";
	}

	file << "]\nk_protocol =  [";
	for(i=0;i<NUM_SEEDS;i++){
		file << "[";
		for(j=0;j<sim_params.total_steps;j++) file << sim_params.k_best_fixed_tau[i*sim_params.total_steps + j] << ", ";
		file << "],";
	}


	file << "]\nj_protocol_secondary =  [";
	for(j=0;j<2*NUMBER_OF_BANGS;j++) file << sim_params.j_best_secondary[j] << ", ";

	file << "]\nk_protocol_secondary =  [";
	for(j=0;j<2*NUMBER_OF_BANGS;j++) file << sim_params.k_best_secondary[j] << ", ";


	file << "]\nbest_evolved_state = [";
	for(i=0;i<2*sim_params.N;i++) file << sim_params.best_evolved_state[i] << ", ";

	file << "]\nbest_mc_result_fixed_tau =   [";
	for(i=0;i<NUM_SEEDS;i++) file << sim_params.best_mc_result_fixed_tau[i] << ", ";

	file << "]\nevolved_state_fixed_tau =   [";
	for(i=0;i<NUM_SEEDS;i++){
		file << "[";
		for(j=0;j<sim_params.N*2;j++) file << sim_params.evolved_state_fixed_tau[i*sim_params.N*2 + j] << ", ";
		file << "],";
	}
	file << "]\n";
	file.close();
}
