#include <fstream>
#include <string>
#include <iostream>

#include "export_data.h"
#include "parameters.h"



void save_mcbf_data_fixed_tau(Simulation_Parameters& sim_params){
	int i,j,seed;
	std::ofstream file;
	std::string type = "MCBF";
	std::string path = make_path(sim_params, type);


	if(sim_params.tau == TAU_INIT_MCBB){
		file.open(path);
		file << "START_PARAMETERS\n";
		file << "DIAG =                       " <<  std::boolalpha << DIAG << "\n";
		file << "NUM_SEEDS =                  " <<  NUM_SEEDS << "\n";
		file << "DISTANCE_LIMIT_MCBF =        " <<  DISTANCE_LIMIT_MCBF << "\n";
		file << "TAU_INIT_MCBF =              " <<  TAU_INIT_MCBF << "\n";
		file << "MAX_TAU_MCBF =               " <<  MAX_TAU_MCBF << "\n";
		file << "TAU_SCALAR_MC =              " <<  TAU_SCALAR_MC << "\n";
		file << "MAX_CHANGE_MCBF_INIT =       " <<  MAX_CHANGE_MCBF_INIT << "\n";
		file << "MIN_CHANGE_MCBF_INIT =       " <<  MIN_CHANGE_MCBF_INIT << "\n";
		file << "ACCEPTANCE_PROB_MC =         " <<  ACCEPTANCE_PROB_MC << "\n";
		file << "TEMP_EXP_DECAY_MC =          " <<  TEMP_EXP_DECAY_MC << "\n";
		file << "BINARY_SEARCH_TAU_LIMIT_MCBF = " <<  BINARY_SEARCH_TAU_LIMIT_MCBF << "\n";
		file << "RANDOM_STATES_MC =           " <<  RANDOM_STATES_MC  << "\n";
		file << "SWEEPS_MC =                  " <<  SWEEPS_MC << "\n";
		file << "TOTAL_STEPS_INIT_MC =        " <<  TOTAL_STEPS_INIT_MC << "\n";
		file << "TEMP_DECAY_ITERATIONS_MC =   " <<  TEMP_DECAY_ITERATIONS_MC  << "\n";
		file << "TEMP_DECAY_LIMIT_MC =        " <<  TEMP_DECAY_LIMIT_MC << "\n";
		file << "MAX_EVOLVE_STEPS_MC =        " <<  MAX_EVOLVE_STEPS_MC << "\n";
		file << "MAX_TAU_STEPS_MCBF =           " <<  MAX_TAU_STEPS_MCBF << "\n";
		file << "ARRAY_SCALAR =               " <<  ARRAY_SCALAR << "\n";
		file.close();
		save_hamiltonian_parameters(sim_params, path);
	}
	file.open(path, std::ios::app);

	file << "\n\n\n\n##############################################################################\n";
	file << "##############################################################################\n";
	file << "tau =               " << sim_params.tau   << "\n";
	file << "total_steps =       " << sim_params.total_steps << "\n";
	file << "time_step =       " << sim_params.time_step << "\n";
	file << "best_E =          " << sim_params.best_E << "\n";
	file << "distance =          " << sim_params.new_distance << "\n";

	file << "\nj_protocol =  [";
	for(i=0;i<NUM_SEEDS;i++){
		file << "[";
		for(j=0;j<sim_params.total_steps;j++) file << sim_params.j_best_fixed_tau[NUMBER_OF_SITES*2*(i*sim_params.total_steps + j)] << ", ";
		file << "],";
	}


	file << "\nk_protocol =  [";
	for(i=0;i<NUM_SEEDS;i++){
		file << "[";
		for(j=0;j<sim_params.total_steps;j++) file << sim_params.k_best_fixed_tau[NUMBER_OF_SITES*2*(i*sim_params.total_steps + j)] << ", ";
		file << "],";
	}

	file << "\nb_protocol =  [";
	for(i=0;i<NUM_SEEDS;i++){
		file << "[";
		for(j=0;j<sim_params.total_steps;j++) file << sim_params.b_best_fixed_tau[NUMBER_OF_SITES*(i*sim_params.total_steps + j)] << ", ";
		file << "],";
	}

	file << "]\nbest_E_fixed_tau =   [";
	for(i=0;i<NUM_SEEDS;i++) file << sim_params.E_array_fixed_tau[i] << ", ";
	file <<"]\n\n";
	file.close();
}

void save_mcbf_data(Simulation_Parameters& sim_params){
	int i,j, cummulative_steps=0;
	std::ofstream file;
	std::string type = "MCBF";
	std::string path = make_path(sim_params, type);
	file.open(path, std::ios::app);

	file << "\n\n\n\n##############################################################################\n";
	file << "##############################################################################\n";
	file << "ground_E =          " << sim_params.ground_E << "\n";
	file << "initial_E =          " << sim_params.initial_E << "\n\n";
	file << "tau_array = [";
	for(i=0;i<sim_params.index;i++) file << sim_params.tau_array[i] << ", ";
	file << "]\nbest_E_array = [";
	for(i=0;i<sim_params.index;i++) file << sim_params.best_E_array[i] << ", ";
	file << "]\ntotal_steps_array = [";
	for(i=0;i<sim_params.index;i++) file << sim_params.total_steps_array[i] << ", ";

	file << "]\nj_best =  [";
	for(i=0;i<sim_params.index;i++){
		file << "[";
		for(j=0;j<sim_params.total_steps_array[i];j++) file << sim_params.j_best[2*NUMBER_OF_SITES*(cummulative_steps+j)] << ", ";
		file << "],";
		cummulative_steps += sim_params.total_steps_array[i];
	}

	cummulative_steps = 0;
	file << "]\nk_best =  [";
	for(i=0;i<sim_params.index;i++){
		file << "[";
		for(j=0;j<sim_params.total_steps_array[i];j++) file << sim_params.k_best[2*NUMBER_OF_SITES*(cummulative_steps+j)] << ", ";
		file << "],";
		cummulative_steps += sim_params.total_steps_array[i];
	}

	cummulative_steps = 0;
	file << "]\nb_best =  [";
	for(i=0;i<sim_params.index;i++){
		file << "[";
		for(j=0;j<sim_params.total_steps_array[i];j++) file << sim_params.b_best[NUMBER_OF_SITES*(cummulative_steps+j)] << ", ";
		file << "],";
		cummulative_steps += sim_params.total_steps_array[i];
	}
	file << "]\n\n";
	file.close();
}



void save_mcbb_data_fixed_tau(Simulation_Parameters& sim_params){
	int i, j;
	std::ofstream file;
	std::string type = "MCBB";
	std::string path = make_path(sim_params, type);

	if(sim_params.tau == TAU_INIT_MCBB){
		file.open(path);
		file << "START_PARAMETERS\n";
		file << "DIAG =                         " <<  std::boolalpha << DIAG << "\n";
		file << "NUM_SEEDS =                    " <<  NUM_SEEDS << "\n";
		file << "DISTANCE_LIMIT_MCBB =          " <<  DISTANCE_LIMIT_MCBB << "\n";
		file << "TAU_INIT_MCBB =                " <<  TAU_INIT_MCBB << "\n";
		file << "MAX_TAU_MCBB =                 " <<  MAX_TAU_MCBB << "\n";
		file << "TAU_SCALAR_MCBB =              " <<  TAU_SCALAR_MCBB << "\n";
		file << "MAX_CHANGE_FRACTION_MCBB =     " <<  MAX_CHANGE_FRACTION_MCBB << "\n";
		file << "MIN_CHANGE_FRACTION_MCBB =     " <<  MIN_CHANGE_FRACTION_MCBB << "\n";
		file << "ACCEPTANCE_PROB_MCBB =         " <<  ACCEPTANCE_PROB_MCBB << "\n";
		file << "TEMP_EXP_DECAY_MCBB =          " <<  TEMP_EXP_DECAY_MCBB << "\n";
		file << "BINARY_SEARCH_TAU_LIMIT_MCBB = " <<  BINARY_SEARCH_TAU_LIMIT_MCBB << "\n";
		file << "RANDOM_STATES_MCBB =           " <<  RANDOM_STATES_MCBB  << "\n";
		file << "NUMBER_OF_BANGS =              " <<  NUMBER_OF_BANGS << "\n";
		file << "SWEEPS_MCBB =                  " <<  SWEEPS_MCBB << "\n";
		file << "TEMP_DECAY_ITERATIONS_MCBB =   " <<  TEMP_DECAY_ITERATIONS_MCBB  << "\n";
		file << "TEMP_DECAY_LIMIT_MCBB =        " <<  TEMP_DECAY_LIMIT_MCBB << "\n";
		file << "MAX_TAU_STEPS_MCBB =           " <<  MAX_TAU_STEPS_MCBB << "\n";
		file.close();
		save_hamiltonian_parameters(sim_params, path);
	}
	file.open(path, std::ios::app);

	file << "\n\n\n\n##############################################################################\n";
	file << "##############################################################################\n";
	file << "tau =               " << sim_params.tau   << "\n";
	file << "best_E =          " << sim_params.best_E << "\n";
	file << "distance =          " << sim_params.new_distance << "\n";

	file << "\nj_protocol =  [";
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

	file << "]\nb_protocol = [";
	for(i=0;i<NUM_SEEDS;i++){
		file << "[";
		for(j=0;j<2*NUMBER_OF_BANGS;j++) file << sim_params.b_best_fixed_tau[i*2*NUMBER_OF_BANGS + j] << ", ";
		file << "],";
	}

	file << "]\nj/k/b =   [";
	for(i=1;i<2*NUMBER_OF_BANGS+1;i++) file << i%2 << ", ";

	file << "]\nbest_E_fixed_tau =   [";
	for(i=0;i<NUM_SEEDS;i++) file << sim_params.E_array_fixed_tau[i] << ", ";
	file <<"]\n\n";
	file.close();
}

void save_mcbb_data(Simulation_Parameters& sim_params){
	int i,j;
	std::ofstream file;
	std::string type = "MCBB";
	std::string path = make_path(sim_params, type);
	file.open(path, std::ios::app);

	file << "\n\n\n\n##############################################################################\n";
	file << "##############################################################################\n";
	file << "ground_E =          " << sim_params.ground_E << "\n";
	file << "initial_E =          " << sim_params.initial_E << "\n\n";
	file << "tau_array = [";
	for(i=0;i<sim_params.index;i++) file << sim_params.tau_array[i] << ", ";
	file << "]\nbest_E_array = [";
	for(i=0;i<sim_params.index;i++) file << sim_params.best_E_array[i] << ", ";

	file << "]\nj_best =  [";
	for(i=0;i<sim_params.index;i++){
		file << "[";
		for(j=0;j<2*NUMBER_OF_BANGS;j++) file << sim_params.j_best[i*2*NUMBER_OF_BANGS + j] << ", ";
		file << "],";
	}
	file << "]\nk_best =  [";
	for(i=0;i<sim_params.index;i++){
		file << "[";
		for(j=0;j<2*NUMBER_OF_BANGS;j++) file << sim_params.k_best[i*2*NUMBER_OF_BANGS + j] << ", ";
		file << "],";
	}
	file << "]\nb_best =  [";
	for(i=0;i<sim_params.index;i++){
		file << "[";
		for(j=0;j<2*NUMBER_OF_BANGS;j++) file << sim_params.b_best[i*2*NUMBER_OF_BANGS + j] << ", ";
		file << "],";
	}	
	file << "]\n\n";
	file.close();
}

void save_adiabatic_data(Simulation_Parameters& sim_params){
	int i;
	std::ofstream file;
	std::string type = "ADIA";
	std::string path = make_path(sim_params, type);

	file.open(path);
	file << "START_PARAMETERS\n";
	file << "DIAG =                       " <<  std::boolalpha << DIAG << "\n";
	file << "DISTANCE_LIMIT_ADIA =        " <<  DISTANCE_LIMIT_ADIA << "\n";
	file << "TAU_INIT_ADIA =              " <<  TAU_INIT_ADIA << "\n";
	file << "MAX_TAU_ADIA =               " <<  MAX_TAU_ADIA << "\n";
	file << "TAU_SCALAR_ADIA =            " <<  TAU_SCALAR_ADIA << "\n";
	file << "TIME_STEP_ADIA =             " <<  TIME_STEP_ADIA << "\n";
	file.close();
	save_hamiltonian_parameters(sim_params, path);
	file.open(path, std::ios::app);

	file << "\n\n\n\n##############################################################################\n";
	file << "##############################################################################\n";
	file << "tau =               " << sim_params.tau   << "\n";
	file << "total_steps =       " << sim_params.total_steps << "\n";
	file << "ground_E =          " << sim_params.ground_E << "\n\n";

	file << "tau_array = [";
	for(i=0;i<sim_params.index;i++) file << sim_params.tau_array[i] << ", ";
	file << "]\nbest_E_array = [";
	for(i=0;i<sim_params.index;i++) file << sim_params.best_E_array[i] << ", ";
	file << "]\n\n";

	file.close();
}



std::string make_path(Simulation_Parameters sim_params, std::string type){
	std::string uni, pbc, gi, gt, fi, ft;

	gi = std::to_string(sim_params.g_initial);
	gi.erase(gi.find_last_not_of('0') + 2, std::string::npos);
	fi = std::to_string(sim_params.f_initial);
	fi.erase(fi.find_last_not_of('0') + 2, std::string::npos);
	gt = std::to_string(sim_params.g_target);
	gt.erase(gt.find_last_not_of('0') + 2, std::string::npos);
	ft = std::to_string(sim_params.f_target);
	ft.erase(ft.find_last_not_of('0') + 2, std::string::npos);

	if (PERIODIC) pbc = 't';
	else pbc = 'f';
	if (UNIFORM_SITES) uni = 't';
	else uni = 'f';

	std::string dir = "../data/" + std::to_string(NX) + "x" +std::to_string(NY) + "/" + std::to_string(sim_params.num_occupants) ;
	std::string file_name = "_occupants/" + type + "___PBC="+pbc+"_UNI="+uni+"_DD="+std::to_string(DEVICE_DIMENSION)+"___gi=" + gi + "_fi=" + fi +"_gt=" + gt +  "_ft="+ ft +".txt";
	std::string path = dir+file_name;

	return path;
}



void save_hamiltonian_parameters(Simulation_Parameters sim_params,std::string path){
	std::ofstream file;
	file.open(path, std::ios::app);

	file << "g_initial =        " << sim_params.g_initial << "\n";
	file << "f_initial =        " << sim_params.f_initial << "\n";
	file << "g_targe t =        " << sim_params.g_target << "\n";
	file << "f_target =         " << sim_params.f_target << "\n";
	file << "j_initial =        " << sim_params.j_initial << "\n";
	file << "k_initial =        " << sim_params.k_initial << "\n";
	file << "b_initial =        " << sim_params.b_initial << "\n";
	file << "j_target =         " << sim_params.j_target << "\n";
	file << "k_target =         " << sim_params.k_target << "\n";
	file << "b_target =         " << sim_params.b_target << "\n";
	file << "PERIODIC =         " << std::boolalpha << PERIODIC << "\n";
	file << "UNIFORM_SITES =    " << std::boolalpha << UNIFORM_SITES << "\n";
	file << "DEVICE_DIMENSION = " << DEVICE_DIMENSION << "\n";
	file << "MAX_PARAM =        " << MAX_PARAM << "\n";
	file << "MIN_PARAM =        " << MIN_PARAM << "\n";
	file << "END_PARAMETERS\n";
	file.close();
}






void save_gap_data(Simulation_Parameters& sim_params, double *f, double *g, double *gap, int L){
	int i;
	std::ofstream file;
	std::string path = "../data/gap.txt";
	file.open(path);

	file << "f = [";
	for(i=0;i<L*L;i++) file << f[i] << ", ";
	file << "]\ng = [";
	for(i=0;i<L*L;i++) file << g[i] << ", ";
	file << "]\ngap = [";
	for(i=0;i<L*L;i++) file << gap[i] << ", ";
	file << "]\n";
	file.close();
}
