#include <string>
#include <fstream>
#include <iostream>



// g++ -o make_submits make_submits.cpp -std=gnu++11



main(int argc, char *argv[]){

	int ji_i, ki_i, jt_i, kt_i, occ_i;
	int num_occupants[1] = {2};
	double j_init[2] = {0.05, 0.95};
	double k_init[2] = {0.05, 0.95};
	double j_targ[7] = {0.05, 0.2, 0.35, 0.5, 0.65, 0.8, 0.95};
	double k_targ[7] = {0.05, 0.2, 0.35, 0.5, 0.65, 0.8, 0.95};
	std::string ji, ki, jt, kt, num, identifier;
	std::ofstream submit_file;

	submit_file.open("submit_file.txt");



	for(occ_i=0;occ_i<(sizeof(num_occupants)/sizeof(int));occ_i++){
		for(ji_i=0;ji_i<(sizeof(j_init)/sizeof(double));ji_i++){
			for(ki_i=0;ki_i<(sizeof(k_init)/sizeof(double));ki_i++){
				for(jt_i=0;jt_i<(sizeof(j_targ)/sizeof(double));jt_i++){
					for(kt_i=0;kt_i<(sizeof(k_targ)/sizeof(double));kt_i++){

						num = std::to_string(num_occupants[occ_i]);
						ji = std::to_string(j_init[ji_i]).substr(0,4);
						ki = std::to_string(k_init[ki_i]).substr(0,4);
						jt = std::to_string(j_targ[jt_i]).substr(0,4);
						kt = std::to_string(k_targ[kt_i]).substr(0,4);


						identifier = "_num_"+num+"__ji_"+ji+"__ki_"+ki+"__jt_"+jt+"__kt_"+kt;
						submit_file << "Executable        = main_4x4\n";
						submit_file << "Arguments         = " << num << " " << ji << " " << ki << " " << jt << " " << kt << "\n";
						submit_file << "Log               = cluster/logs/_" << identifier  << ".log\n";
						submit_file << "Output            = cluster/outputs/_" << identifier << "\n";
						submit_file << "request_cpus      = 4\n";
						//submit_file << "request_memory    = 20 GB\n";
						submit_file << "queue\n";
	}}}}}
}
