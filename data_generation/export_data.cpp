#include <fstream>
#include <string>
#include <iostream>

#include "export_data.h"
#include "parameters.h"



void save_mcbf_data(Simulation_Parameters& sim_params){
	int i;
	std::ofstream file;
	std::string type = "MCBF";
	std::string path = make_path(sim_params, type);


	if(sim_params.seed == 1)	{
		file.open(path);
		file << "START_PARAMETERS\n";
		file << "DIAG =                       " <<  std::boolalpha << DIAG << "\n";
		file << "SEED_TOTAL =                 " <<  SEED_TOTAL << "\n";
		file << "DIFFERENCE_LIMIT_MC =        " <<  DIFFERENCE_LIMIT_MC << "\n";
		file << "TAU_INIT_MC =                " <<  TAU_INIT_MC << "\n";
		file << "MAX_TAU_MC =                 " <<  MAX_TAU_MC << "\n";
		file << "TAU_SCALAR_MC =              " <<  TAU_SCALAR_MC << "\n";
		file << "MAX_CHANGE_MC_INIT =         " <<  MAX_CHANGE_MC_INIT << "\n";
		file << "ACCEPTANCE_PROB_MC =         " <<  ACCEPTANCE_PROB_MC << "\n";
		file << "TEMP_EXP_DECAY_MC =          " <<  TEMP_EXP_DECAY_MC << "\n";
		file << "BINARY_SEARCH_TAU_LIMIT_MC = " <<  BINARY_SEARCH_TAU_LIMIT_MC << "\n";
		file << "RANDOM_STATES_MC =           " <<  RANDOM_STATES_MC  << "\n";
		file << "SWEEPS_MC =                  " <<  SWEEPS_MC << "\n";
		file << "TOTAL_STEPS_INIT_MC =        " <<  TOTAL_STEPS_INIT_MC << "\n";
		file << "TEMP_DECAY_ITERATIONS_MC =   " <<  TEMP_DECAY_ITERATIONS_MC  << "\n";
		file << "TEMP_DECAY_LIMIT_MC =        " <<  TEMP_DECAY_LIMIT_MC << "\n";
		file << "MAX_EVOLVE_STEPS_MC =        " <<  MAX_EVOLVE_STEPS_MC << "\n";
		file << "MAX_TAU_STEPS_MC =           " <<  MAX_TAU_STEPS_MC << "\n";
		file << "ARRAY_SCALAR =               " <<  ARRAY_SCALAR << "\n";
		file.close();
		save_hamiltonian_parameters(sim_params, path);
	}
	file.open(path, std::ios::app);

	file << "\n\n\n\n##############################################################################\n";
	file << "##############################################################################\n";
	file << "seed =              " << sim_params.seed << "\n";
	file << "tau =               " << sim_params.tau   << "\n";
	file << "total_steps =       " << sim_params.total_steps << "\n";
	file << "ground_E =          " << sim_params.ground_E << "\n\n";

	file << "j_protocol = [";
	for(i=0;i<sim_params.total_steps_cummulative_sum[sim_params.index-1];i++) file << sim_params.j_best[i*NUMBER_OF_SITES*2] << ", ";
	file <<"]\nk_protocol = [";
	for(i=0;i<sim_params.total_steps_cummulative_sum[sim_params.index-1];i++) file << sim_params.k_best[i*NUMBER_OF_SITES*2] << ", ";
	file <<"]\nb_protocol = [";
	for(i=0;i<sim_params.total_steps_cummulative_sum[sim_params.index-1];i++) file << sim_params.b_best[i*NUMBER_OF_SITES] << ", ";
	file <<"]\n\n";

	file <<  "tau_array =  [";
	for(i=0;i< sim_params.index;i++) file << sim_params.tau_array[i] << ", ";
	file << "]\nbest_E_array =  [";
	for(i=0;i< sim_params.index;i++) file << sim_params.best_E_array[i] << ", ";
	file << "]\ntotal_steps_cummulative_sum = [";
	for(i=0;i< sim_params.index;i++) file << sim_params.total_steps_cummulative_sum[i] << ", ";
	file << "]\n\n";

	file.close();
}



void save_mcbb_data(Simulation_Parameters& sim_params){
	int i;
	std::ofstream file;
	std::string type = "MCBB";
	std::string path = make_path(sim_params, type);

	if(sim_params.seed == 1){
		file.open(path);
		file << "START_PARAMETERS\n";
		file << "DIAG =                         " <<  std::boolalpha << DIAG << "\n";
		file << "SEED_TOTAL =                   " <<  SEED_TOTAL << "\n";
		file << "DIFFERENCE_LIMIT_MCBB =        " <<  DIFFERENCE_LIMIT_MCBB << "\n";
		file << "TAU_INIT_MCBB =                " <<  TAU_INIT_MCBB << "\n";
		file << "MAX_TAU_MCBB =                 " <<  MAX_TAU_MCBB << "\n";
		file << "TAU_SCALAR_MCBB =              " <<  TAU_SCALAR_MCBB << "\n";
		file << "CHANGE_FRACTION_MCBB =         " <<  CHANGE_FRACTION_MCBB << "\n";
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
	file << "seed =              " << sim_params.seed << "\n";
	file << "tau =               " << sim_params.tau   << "\n";
	file << "ground_E =          " << sim_params.ground_E << "\n\n";

	file << "\nj_protocol =  [";
	for(i=0;i<2*NUMBER_OF_BANGS*sim_params.index;i++) file << sim_params.j_best[i] << ", ";
	file << "]\nk_protocol = [";
	for(i=0;i<2*NUMBER_OF_BANGS*sim_params.index;i++) file << sim_params.k_best[i] << ", ";
	file << "]\nb_protocol = [";
	for(i=0;i<2*NUMBER_OF_BANGS*sim_params.index;i++) file << sim_params.b_best[i] << ", ";
	file << "]\nj/k/b =   [";
	for(i=1;i<2*NUMBER_OF_BANGS+1;i++) file << i%2 << ", ";
	file <<"]\n\n";


	file << "tau_array = [";
	for(i=0;i<sim_params.index;i++) file << sim_params.tau_array[i] << ", ";
	file << "]\nbest_E_array = [";
	for(i=0;i<sim_params.index;i++) file << sim_params.best_E_array[i] << ", ";
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
	file << "DIFFERENCE_LIMIT_ADIA =      " <<  DIFFERENCE_LIMIT_ADIA << "\n";
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
