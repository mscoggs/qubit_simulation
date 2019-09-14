#include <string>
#include <fstream>
#include <iostream>



// g++ -o make_submits make_submits.cpp -std=gnu++11



main(int argc, char *argv[]){

	int gi_i, fi_i, gt_i, ft_i, occ_i;
	int num_occupants[4] = {2,7,3,4};
	double g_init[2] = {0.1, 1};
	double f_init[2] = {0.1, 1};
	double g_targ[3] = {0.1, 0.5, 1};
	double f_targ[3] = {0.1, 0.5, 1};
	std::string submit_file_name, gi, fi, gt, ft, num, identifier;
	std::ofstream submit_file, run_file;
	run_file.open("../run_jobs.sh");
	run_file << "#!/bin/bash\n\n";

	run_file << "# chmod u+x run_jobs.sh\n";
	run_file << "# ./run_jobs.sh\n\n\n";
	for(occ_i=0;occ_i<(sizeof(num_occupants)/sizeof(int));occ_i++){
		for(gi_i=0;gi_i<(sizeof(g_init)/sizeof(double));gi_i++){
			for(fi_i=0;fi_i<(sizeof(f_init)/sizeof(double));fi_i++){
				for(gt_i=0;gt_i<(sizeof(g_targ)/sizeof(double));gt_i++){
					for(ft_i=0;ft_i<(sizeof(f_targ)/sizeof(double));ft_i++){
					
						num = std::to_string(num_occupants[occ_i]);
						gi = std::to_string(g_init[gi_i]).substr(0,3);
						fi = std::to_string(f_init[fi_i]).substr(0,3);
						gt = std::to_string(g_targ[gt_i]).substr(0,3);
						ft = std::to_string(f_targ[ft_i]).substr(0,3);

						identifier = "_num_"+num+"__gi_"+gi+"__fi_"+fi+"__gt_"+gt+"__ft_"+ft;
						submit_file_name = "submits/submit_" + identifier  +".job";
						
						submit_file.open(submit_file_name);
						submit_file << "Executable        = main\n";
						submit_file << "Arguments         = " << num << " " << gi << " " << fi << " " << gt << " " << ft << "\n";
						submit_file << "Log               = cluster/logs/_" << identifier  << ".log\n";
						submit_file << "Output            = cluster/outputs/_" << identifier << "\n";
						submit_file << "request_cpus      = 2\n";
						submit_file << "Queue\n";
						submit_file.close();

						run_file << "condor_submit cluster/" << submit_file_name << "\n";
	}}}}}
}
