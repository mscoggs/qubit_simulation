#include <string>
#include <fstream>
#include <iostream>



// g++ -o make_submits make_submits.cpp -std=gnu++11



main(int argc, char *argv[]){

	int ji_i, fi_i, jt_i, kt_i, occ_i;
	int num_occupants[2] = {3,4};
	double j_init[2] = {0.1, 1};
	double k_init[2] = {0.1, 1};
	double j_targ[3] = {0.1, 0.5, 1};
	double k_targ[3] = {0.1, 0.5, 1};
	std::string ji, fi, jt, kt, num, identifier;
	std::ofstream submit_file;


	submit_file.open("../submit_file.txt");



	for(occ_i=0;occ_i<(sizeof(num_occupants)/sizeof(int));occ_i++){
		for(ji_i=0;ji_i<(sizeof(g_init)/sizeof(double));ji_i++){
			for(ki_i=0;ki_i<(sizeof(f_init)/sizeof(double));ki_i++){
				for(jt_i=0;jt_i<(sizeof(g_targ)/sizeof(double));jt_i++){
					for(kt_i=0;kt_i<(sizeof(f_targ)/sizeof(double));kt_i++){

						num = std::to_string(num_occupants[occ_i]);
						ji = std::to_string(g_init[ji_i]).substr(0,3);
						ki = std::to_string(f_init[ki_i]).substr(0,3);
						jt = std::to_string(g_targ[jt_i]).substr(0,3);
						kt = std::to_string(f_targ[kt_i]).substr(0,3);


						identifier = "_num_"+num+"__ji_"+ji+"__ki_"+ki+"__jt_"+jt+"__kt_"+kt;

						submit_file << "Executable        = main\n";
						submit_file << "Arguments         = " << num << " " << ji << " " << ki << " " << jt << " " << kt << "\n";
						submit_file << "Log               = cluster/logs/_" << identifier  << ".log\n";
						submit_file << "Output            = cluster/outputs/_" << identifier << "\n";
						submit_file << "request_cpus      = 2\n";
						//submit_file << "request_memory    = 20 GB\n";
						submit_file << "queue\n";
	}}}}}
}
