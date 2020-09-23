#include <fstream>
#include <string>
#include <iostream>

#include "write_data.h"
#include "parameters.h"



void save_mcbf_data(Simulation_Parameters& sim_params){
	int i,j,seed;
	std::ofstream file;
	std::string type = "MCBF";
	std::string path = make_path(sim_params, type);


	if(sim_params.tau == TAU_INIT || sim_params.tau == 0.0){
		file.open(path);
		file << "START_PARAMETERS\n";
		file << "MAX_CHANGE_MCBF =       " <<  MAX_CHANGE_MCBF << "\n";
		file << "MIN_CHANGE_MCBF =       " <<  MIN_CHANGE_MCBF << "\n";
		file << "SWEEPS_MCBF =                  " <<  SWEEPS_MCBF << "\n";
		file << "TOTAL_STEPS_INIT_MCBF =        " <<  TOTAL_STEPS_INIT_MCBF << "\n";
		file << "MAX_EVOLVE_STEPS_MCBF =        " <<  MAX_EVOLVE_STEPS_MCBF << "\n";
		file << "ARRAY_SCALAR =               " <<  ARRAY_SCALAR << "\n";
		file.close();
		save_hamiltonian_parameters(sim_params, path);
	}
	file.open(path, std::ios::app);

	file << "##############################################################################\n";
	file << "clock_duration  =         " << sim_params.duration  << "\n";
	file << "tau =               " << sim_params.tau   << "\n";
	file << "total_steps =       " << sim_params.total_steps << "\n";
	file << "time_step =       " << sim_params.time_step << "\n";
	file << "best_mc_result =          " << sim_params.best_mc_result << "\n";
	file << "evolved_target_dot_squared =   " << sim_params.evolved_target_dot_squared << "\n";
	file << "best_evolved_state = [";
	for(i=0;i<2*sim_params.N;i++) file << sim_params.best_evolved_state[i] << ", ";
	file << "]\n";
	file << "distance =          " << sim_params.new_distance << "\n";

	file << "j_protocol =  [";
	for(i=0;i<NUM_SEEDS;i++){
		file << "[";
		for(j=0;j<sim_params.total_steps;j++) file << sim_params.j_best_fixed_tau[NUMBER_OF_SITES*2*(i*sim_params.total_steps + j)] << ", ";
		file << "],";
	}

	file << "]\nk_protocol =  [";
	for(i=0;i<NUM_SEEDS;i++){
		file << "[";
		for(j=0;j<sim_params.total_steps;j++) file << sim_params.k_best_fixed_tau[NUMBER_OF_SITES*2*(i*sim_params.total_steps + j)] << ", ";
		file << "],";
	}

	// file << "]\nb_protocol =  [";
	// for(i=0;i<NUM_SEEDS;i++){
	// 	file << "[";
	// 	for(j=0;j<sim_params.total_steps;j++) file << sim_params.b_best_fixed_tau[NUMBER_OF_SITES*(i*sim_params.total_steps + j)] << ", ";
	// 	file << "],";
	// }

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




void save_mcbb_data(Simulation_Parameters& sim_params){
	int i, j;
	std::ofstream file;
	std::string type = "MCBB";
	std::string path = make_path(sim_params, type);

	if(sim_params.tau == TAU_INIT || sim_params.tau == 0.0){
		file.open(path);
		file << "START_PARAMETERS\n";
		file << "MAX_CHANGE_FRACTION_MCBB =     " <<  MAX_CHANGE_FRACTION_MCBB << "\n";
		file << "MIN_CHANGE_FRACTION_MCBB =     " <<  MIN_CHANGE_FRACTION_MCBB << "\n";
		file << "NUMBER_OF_BANGS =              " <<  NUMBER_OF_BANGS << "\n";
		file << "SWEEPS_MCBB =                  " <<  SWEEPS_MCBB << "\n";
		file.close();
		save_hamiltonian_parameters(sim_params, path);
	}
	file.open(path, std::ios::app);

	file << "##############################################################################\n";
	file << "clock_duration  =         " << sim_params.duration  << "\n";
	file << "tau =               " << sim_params.tau   << "\n";
	file << "best_mc_result =          " << sim_params.best_mc_result << "\n";
	file << "evolved_target_dot_squared =   " << sim_params.evolved_target_dot_squared << "\n";
	file << "best_evolved_state = [";
	for(i=0;i<2*sim_params.N;i++) file << sim_params.best_evolved_state[i] << ", ";
	file << "]\n";
	file << "distance =          " << sim_params.new_distance << "\n";

	file << "j_protocol =  [";
	for(i=0;i<NUM_SEEDS;i++){
		file << "[";
		for(j=0;j<2*NUMBER_OF_BANGS;j++) file << sim_params.j_best_fixed_tau[i*2*NUMBER_OF_BANGS + j] << ", ";
		file << "],";
	}

	file << "]\nk_protocol = [";
	for(i=0;i<NUM_SEEDS;i++){
		file << "[";
		for(j=0;j<2*NUMBER_OF_BANGS;j++) file << sim_params.k_best_fixed_tau[i*2*NUMBER_OF_BANGS + j] << ", ";
		file << "],";
	}

	// file << "]\nb_protocol = [";
	// for(i=0;i<NUM_SEEDS;i++){
	// 	file << "[";
	// 	for(j=0;j<2*NUMBER_OF_BANGS;j++) file << sim_params.b_best_fixed_tau[i*2*NUMBER_OF_BANGS + j] << ", ";
	// 	file << "],";
	// }

	file << "]\nj/k/b =   [";
	for(i=1;i<2*NUMBER_OF_BANGS+1;i++) file << i%2 << ", ";

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




void save_mcdb_data(Simulation_Parameters& sim_params){
	int i, j;
	std::ofstream file;
	std::string type = "MCDB";
	std::string path = make_path(sim_params, type);

	if(sim_params.tau == TAU_INIT || sim_params.tau == 0.0){
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
		file.close();
		save_hamiltonian_parameters(sim_params, path);
	}
	file.open(path, std::ios::app);

	file << "##############################################################################\n";
	file << "clock_duration  =         " << sim_params.duration  << "\n";
	file << "tau =                     " << sim_params.tau   << "\n";
	file << "total_steps =             " << sim_params.total_steps << "\n";
	file << "time_step =               " << sim_params.time_step << "\n";
	file << "best_mc_result =          " << sim_params.best_mc_result << "\n";
	file << "evolved_target_dot_squared =  " << sim_params.evolved_target_dot_squared << "\n";
	file << "best_evolved_state = [";
	for(i=0;i<2*sim_params.N;i++) file << sim_params.best_evolved_state[i] << ", ";
	file << "]\n";
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

	// file << "]\nb_protocol =  [";
	// for(i=0;i<NUM_SEEDS;i++){
	// 	file << "[";
	// 	for(j=0;j<sim_params.total_steps;j++) file << sim_params.b_best_fixed_tau[i*sim_params.total_steps + j] << ", ";
	// 	file << "],";
	// }

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




void save_adiabatic_data(Simulation_Parameters& sim_params){
	int i;
	std::ofstream file;
	std::string type = "ADIA";
	std::string path = make_path(sim_params, type);

	if(sim_params.tau == TAU_INIT || sim_params.tau == 0.0){
		file.open(path);
		file << "START_PARAMETERS\n";
		file << "MAX_TAU_ADIA =             " <<  MAX_TAU_ADIA << "\n";
		file << "TOTAL_STEPS_ADIA =             " <<  TOTAL_STEPS_ADIA << "\n";
		file.close();
		save_hamiltonian_parameters(sim_params, path);
	}
	file.open(path, std::ios::app);

	file << "##############################################################################\n";
	file << "tau =               " << sim_params.tau   << "\n";
	file << "total_steps =       " << sim_params.total_steps << "\n";
	file << "time_step  =        " << sim_params.time_step  << "\n";
	file << "best_mc_result  =           " << sim_params.best_mc_result  << "\n";
	file << "evolved_target_dot_squared =   " << sim_params.evolved_target_dot_squared << "\n";
	file << "best_evolved_state = [";
	for(i=0;i<2*sim_params.N;i++) file << sim_params.best_evolved_state[i] << ", ";
	file << "]\n";
	file << "distance  =         " << sim_params.new_distance  << "\n";
	file << "clock_duration  =         " << sim_params.duration  << "\n";
	file.close();
}



std::string make_path(Simulation_Parameters sim_params, std::string type){
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
	std::string dir = "../data_binary_big/" + std::to_string(NX) + "x" +std::to_string(NY) + "/" + std::to_string(sim_params.num_occupants) ;
	std::string file_name = "_occupants/" + type + "___PBC="+pbc+"_UNI="+uni+"_DD="+std::to_string(DEVICE_DIMENSION)+"___ri=" + ri +  "_rt="+ rt +".txt";
	std::string path = dir+file_name;

	return path;
}



void save_hamiltonian_parameters(Simulation_Parameters sim_params,std::string path){
	int i;
	std::ofstream file;


	file.open(path, std::ios::app);
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
	file << "DISTANCE_LIMIT =          " <<  DISTANCE_LIMIT << "\n";
	file << "TAU_INIT =                " <<  TAU_INIT << "\n";
	file << "MAX_TAU =                 " <<  MAX_TAU << "\n";
	file << "TAU_SCALAR =              " <<  TAU_SCALAR << "\n";
	file << "TAU_SCALAR_TINY =         " <<  TAU_SCALAR_TINY << "\n";
	file << "TAU_SCALAR_BIG =          " <<  TAU_SCALAR_BIG << "\n";
	file << "ACCEPTANCE_PROB =         " <<  ACCEPTANCE_PROB << "\n";
	file << "TEMP_EXP_DECAY =          " <<  TEMP_EXP_DECAY << "\n";
	file << "MIN_TEMP_FRACTION =       " <<  MIN_TEMP_FRACTION << "\n";
	file << "TEMP_DECAY_ITERATIONS =   " <<  TEMP_DECAY_ITERATIONS  << "\n";
	file << "ZERO_TEMP_ITREATIONS =    " <<  ZERO_TEMP_ITERATIONS  << "\n";
	file << "RANDOM_STATES =           " <<  RANDOM_STATES  << "\n";
	file << "INIT_OVERLAP_LIMIT =      " <<  INIT_OVERLAP_LIMIT  << "\n";
	file << "END_PARAMETERS\n";
	file.close();
}






void save_gap_data(Simulation_Parameters& sim_params, double *j, double *k, double *gap, int L){
	int i;
	std::ofstream file;
	std::string path = "../data/gap.txt";
	file.open(path);

	file << "j = [";
	for(i=0;i<L*L;i++) file << j[i] << ", ";
	file << "]\nk = [";
	for(i=0;i<L*L;i++) file << k[i] << ", ";
	file << "]\ngap = [";
	for(i=0;i<L*L;i++) file << gap[i] << ", ";
	file << "]\n";
	file.close();
}
