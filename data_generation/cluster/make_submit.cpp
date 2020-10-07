#include <string>
#include <math.h>
#include <fstream>
#include <iostream>



// g++ -o make_submits make_submits.cpp -std=gnu++11



main(int argc, char *argv[]){

	int ji_i, ki_i, jt_i, kt_i, occ_i;
	int num_occupants[1] = {3};
	std::string ji, ki, jt, kt, num, identifier;
	double j_init, k_init, j_targ, k_targ;
	double log_ri_min = -3, log_ri_max = 3;
	double log_rt_min = -3, log_rt_max = 3;
	int ri_points = 100;
	int rt_points = 100;
	double log_ri, log_rt;
	std::ofstream submit_file;

	submit_file.open("submit_file.txt");



	for(occ_i=0;occ_i<(sizeof(num_occupants)/sizeof(int));occ_i++){
		for(log_ri = log_ri_min;log_ri<=log_ri_max;log_ri+=(log_ri_max-log_ri_min)/ri_points){
			for(log_rt = log_rt_min;log_rt<=log_rt_max;log_rt+=(log_rt_max-log_rt_min)/rt_points){
				double ri = exp(log_ri), rt = exp(log_rt);
			       	k_init = 0.5;
				j_init = ri*k_init;
				while(j_init >= 1){
					j_init = j_init/2.0;
					k_init = k_init/2.0;
				}
			       	k_targ = 0.5;
				j_targ = rt*k_targ;
				while(j_targ >= 1){
					j_targ = j_targ/2.0;
					k_targ = k_targ/2.0;
				}

				num = std::to_string(num_occupants[occ_i]);
				ji = std::to_string(j_init).substr(0,5);
				ki = std::to_string(k_init).substr(0,5);
				jt = std::to_string(j_targ).substr(0,5);
				kt = std::to_string(k_targ).substr(0,5);

				identifier = "_num_"+num+"__ji_"+ji+"__ki_"+ki+"__jt_"+jt+"__kt_"+kt;
				submit_file << "Executable        = main_3x3_3\n";
				submit_file << "Arguments         = " << num << " " << ji << " " << ki << " " << jt << " " << kt << "\n";
				//submit_file << "Log               = cluster/logs/_" << identifier  << ".log\n";
				//submit_file << "Output            = cluster/outputs/_" << identifier << "\n";
	//						submit_file << "request_cpus      = 4\n";
				//submit_file << "request_memory    = 20 GB\n";
				submit_file << "queue\n";
	}}}
	}
