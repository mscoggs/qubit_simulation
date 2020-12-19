#include <stdio.h>

#include "parameters.h"
#include "print.h"



void print_state(double* state,int N){
	int i;
	printf("[");
	for(i=0;i<N;i++) printf("%04.3f+%04.3fi; ", state[2*i], state[2*i+1]);
	printf("]\n");
}



void print_hamiltonian_complex(double* hamiltonian, int N){
	int i,j;
	printf("\nPrinting the %ix%i Hamiltonian=\n[", N,N);
	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++) printf("%09.5f+%09.5fi  ",(hamiltonian[2*(j*N+i)]+0.0), hamiltonian[2*(j*N+i)+1]);
		if (i==N-1) printf("]");
		else printf(";\n ");
	}
	printf("\n");
}



void print_hamiltonian_real(double* hamiltonian, int N){
	int i,j;
	printf("\nPrinting the %ix%i Hamiltonian=\n[", N,N);
	for (i=0;i<N;i++){
		for (j=0;j<N;j++) printf("%09.5f ",(hamiltonian[j*N+i])+0.0);
		if (i==N-1) printf("]");
		else printf("\n ");
	}
	printf("\n");
}



void print_arrays_mcbb(double* j_times, double* k_times, double* b_times){
	int i;
	printf("\nj_time=[");

	for (i=0; i<2*NUMBER_OF_BANGS; i++) printf(" %5.4f |", j_times[i]);

	printf("]\nk_times=[");
	for (i=0; i<2*NUMBER_OF_BANGS; i++) printf(" %5.4f |", k_times[i]);

	//printf("]\nb_times=[");
	//for (i=0; i<2*NUMBER_OF_BANGS; i++) printf(" %5.4f |", b_times[i]);
	printf("\n");
}

void print_arrays_mcdb(double* j_times, double* k_times, double* b_times, int size){
	int i;
	printf("\nj_time=[");

	for (i=0; i<size; i++) printf(" %5.4f |", j_times[i]);

	printf("]\nk_times=[");
	for (i=0; i<size; i++) printf(" %5.4f |", k_times[i]);

	//printf("]\nb_times=[");
	//for (i=0; i<size; i++) printf(" %5.4f |", b_times[i]);
	printf("\n");
}


void print_arrays_mcbf(double *j_array, double *k_array, double *b_array, int total_steps){
	int i,j;

	if (UNIFORM_SITES){
		printf("\nj_array= %iX%i (stepsXsites)\n[", total_steps, NUMBER_OF_SITES*2);
		for(i=0;i<total_steps;i++) printf("%5.2f ",j_array[i*NUMBER_OF_SITES*2]);
		printf(";]\n");

		printf("k_array= %iX%i (stepsXsites)\n[", total_steps, NUMBER_OF_SITES*2);
		for(i=0;i<total_steps;i++) printf("%5.2f ",k_array[i*NUMBER_OF_SITES*2]);
		printf(";]\n");

		printf("b_array= %iX%i (stepsXsites)\n[", total_steps, NUMBER_OF_SITES);
		for(i=0;i<total_steps;i++) printf("%5.2f ",b_array[i*NUMBER_OF_SITES]);
		printf(";]\n");
	}

	else{
		printf("\nj_array= %iX%i (stepsXsites)\n[", total_steps, NUMBER_OF_SITES*2);
		for(i=0;i<total_steps;i++){
			for(j=0;j<2*NUMBER_OF_SITES;j++) printf("%5.2f ",j_array[i*NUMBER_OF_SITES*2+j]);
			if (i==total_steps-1) printf("]\n");
			else printf(";\n ");
		}

		printf("k_array= %iX%i (stepsXsites)\n[", total_steps, NUMBER_OF_SITES*2);
		for(i=0;i<total_steps;i++){
			for(j=0;j<2*NUMBER_OF_SITES;j++) printf("%5.2f ",k_array[i*NUMBER_OF_SITES*2+j]);
			if (i==total_steps-1) printf("]\n");
			else printf(";\n ");
		}

		printf("b_array= %iX%i (stepsXsites)\n[", total_steps, NUMBER_OF_SITES);
		for(i=0;i<total_steps;i++){
			for(j=0;j<NUMBER_OF_SITES;j++) printf("%5.2f ",b_array[i*NUMBER_OF_SITES+j]);
			if (i==total_steps-1) printf("]\n");
			else printf(";\n ");
		}
		printf("\n");
	}
}



void print_mc_results(Simulation_Parameters& sim_params){
	printf("############################################################################\n");
	printf("BEST RESULT=    %9.6f\n",sim_params.best_mc_result);
	printf("OLD DISTANCNCE= %9.6f\n",sim_params.old_distance);
	printf("NEW DISTANCNCE= %9.6f\n",sim_params.new_distance);
	printf("############################################################################\n\n\n");
}



void print_adiabatic_info(Simulation_Parameters& sim_params){
		printf("############################################################################\n");
		printf("########################## THE ADIABATIC METHOD ############################\n");
		printf("############################################################################\n");
		printf("|| OCCUPANTS:       %4i || DIMENSION:        %4i ||                     ||\n", sim_params.num_occupants, sim_params.N);
		printf("|| TAU_MAX:       %4.3f || TOTAL_STEPS:      %4i || TIME_STEP    %4.4f ||\n", double(MAX_TAU), sim_params.total_steps, sim_params.time_step);
		printf("|| J_INITIAL:     %4.4f || K_INITIAL:      %4.4f || B_INITIAL:   %4.4f ||\n", sim_params.j_initial,     sim_params.k_initial,   sim_params.b_initial);
		printf("|| J_TARGET:      %4.4f || K_TARGET:       %4.4f || B_TARGET:    %4.4f ||\n", sim_params.j_target,      sim_params.k_target,    sim_params.b_target);
		printf("|| GROUND_E:     %4.4f || INITIAL_E:     %4.4f ||                     ||\n", sim_params.ground_E,      sim_params.initial_E);
		printf("\nINITIAL_STATE: "), print_state(sim_params.init_state, sim_params.N);
		printf("############################################################################\n");
}



void print_mcbf_info(Simulation_Parameters& sim_params){
		printf("############################################################################\n");
		printf("################### THE  MONTE-CARLO BRUTE-FORCE METHOD ####################\n");
		printf("############################################################################\n");
		printf("|| OCCUPANTS:       %4i || DIMENSION:        %4i || SEED:          %4i ||\n", sim_params.num_occupants, sim_params.N, sim_params.seed);
		printf("|| TAU_MAX:       %4.4f || TAU:            %4.4f || TIME_STEP    %4.4f ||\n", double(MAX_TAU), sim_params.tau, sim_params.time_step);
		printf("|| TOTAL_STEPS:     %4i || TOTAL_SWEEPS:     %4i ||                     ||\n", sim_params.total_steps, SWEEPS_MCBF*sim_params.total_steps*sim_params.sweeps_multiplier);
		printf("|| TEMPERATURE:   %4.4f || TEMP_DECAYS:      %4i ||                     ||\n", sim_params.temperature, TEMP_DECAY_ITERATIONS);
		printf("|| J_INITIAL:     %4.4f || K_INITIAL:      %4.4f || B_INITIAL:   %4.4f ||\n", sim_params.j_initial,     sim_params.k_initial,   sim_params.b_initial);
		printf("|| J_TARGET:      %4.4f || K_TARGET:       %4.4f || B_TARGET:    %4.4f ||\n", sim_params.j_target,      sim_params.k_target,    sim_params.b_target);
		printf("|| GROUND_E:     %4.4f || INITIAL_E:     %4.4f ||                     ||\n", sim_params.ground_E,      sim_params.initial_E);
		printf("\nINITIAL_STATE: "), print_state(sim_params.init_state, sim_params.N);
		printf("############################################################################\n");
}



void print_mcbb_info(Simulation_Parameters& sim_params){
		printf("############################################################################\n");
		printf("################### THE  MONTE-CARLO BANG-BANG METHOD ######################\n");
		printf("############################################################################\n");
		printf("|| OCCUPANTS:       %4i || DIMENSION:        %4i || SEED:          %4i ||\n", sim_params.num_occupants, sim_params.N, sim_params.seed);
		printf("|| TAU_MAX:       %4.4f || TAU:            %4.4f || TOTAL_SWEEPS:  %4i ||\n", double(MAX_TAU), sim_params.tau, SWEEPS_MCBB*NUMBER_OF_BANGS*sim_params.sweeps_multiplier);
		printf("|| TEMPERATURE:   %4.4f || TEMP_DECAYS:      %4i ||                     ||\n", sim_params.temperature, TEMP_DECAY_ITERATIONS);
		printf("|| J_INITIAL:     %4.4f || K_INITIAL:      %4.4f || B_INITIAL:   %4.4f ||\n", sim_params.j_initial,     sim_params.k_initial,   sim_params.b_initial);
		printf("|| J_TARGET:      %4.4f || K_TARGET:       %4.4f || B_TARGET:    %4.4f ||\n", sim_params.j_target,      sim_params.k_target,    sim_params.b_target);
		printf("|| GROUND_E:     %4.4f || INITIAL_E:     %4.4f ||                     ||\n", sim_params.ground_E,      sim_params.initial_E);
		printf("\nINITIAL_STATE: "), print_state(sim_params.init_state, sim_params.N);
		printf("############################################################################\n");
}




void print_mcdb_info(Simulation_Parameters& sim_params){
		printf("############################################################################\n");
		printf("################# THE  MONTE-CARLO DISCRETE-BANG METHOD ####################\n");
		printf("############################################################################\n");
		printf("|| OCCUPANTS:       %4i || DIMENSION:        %4i || SEED:          %4i ||\n", sim_params.num_occupants, sim_params.N, sim_params.seed);
		printf("|| TAU_MAX:       %4.4f || TAU:            %4.4f || TIME_STEP    %4.4f ||\n", double(MAX_TAU), sim_params.tau, sim_params.time_step);
		printf("|| TOTAL_STEPS:     %4i || TOTAL_SWEEPS:     %4i ||                     ||\n", sim_params.total_steps, sim_params.total_sweeps*sim_params.sweeps_multiplier);
		printf("|| TEMPERATURE:   %4.4f || TEMP_DECAYS:      %4i ||                     ||\n", sim_params.temperature, TEMP_DECAY_ITERATIONS);
		printf("|| J_INITIAL:     %4.4f || K_INITIAL:      %4.4f || B_INITIAL:   %4.4f ||\n", sim_params.j_initial,     sim_params.k_initial,   sim_params.b_initial);
		printf("|| J_TARGET:      %4.4f || K_TARGET:       %4.4f || B_TARGET:    %4.4f ||\n", sim_params.j_target,      sim_params.k_target,    sim_params.b_target);
		printf("|| TARGET:        %4.4f || INITIAL:        %4.4f || DISTANCE:    %4.4f ||\n", sim_params.ground_E,      sim_params.initial_E, sim_params.new_distance);
//		printf("\nINITIAL_STATE: "), print_state(sim_params.init_state, sim_params.N);
		printf("############################################################################\n");
}
