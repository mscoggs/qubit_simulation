#include <string>
#include <cmath>
#include <fstream>
#include <iostream>



// g++ -o make_q_jobs make_q_jobs.cpp -std=gnu++11



main(int argc, char *argv[]){

	int num_jobs = atoi(argv[1]), job_num;
	double jt=0;
	std::ofstream job_file;
	job_file.open("submit_q_file.txt");
	job_file << "#!/bin/bash\n";

	

	for(job_num=0;job_num<num_jobs;job_num++){
		jt = floor((job_num*1000.0)/num_jobs)/1000.0;
		printf("jt: %f\n", jt);
		job_file << "./main 2 0.8 0.2 0.2 " <<std::to_string(jt) << "&\n";

	}
}
