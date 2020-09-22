#include <gsl/gsl_rng.h>
#include <cstring>
#include <algorithm>    // std::min

#include "check.h"
#include "write_data.h"
#include "hamiltonian.h"
#include "linear_algebra.h"
#include "mcdb.h"
#include "operations.h"
#include "parameters.h"
#include "print.h"




void mcdb_method(Simulation_Parameters& sim_params){
	int i;


	sim_params.init_mcdb_params();

	if(check_commutator(sim_params.N, sim_params.ham_initial, sim_params.ham_target) || (sim_params.init_target_dot_squared > INIT_OVERLAP_LIMIT && USE_ENERGY_DISTANCE) || (1-sim_params.init_target_dot_squared < DISTANCE_LIMIT && !USE_ENERGY_DISTANCE)){
		sim_params.tau = 0.0, sim_params.new_distance = 0.0, sim_params.best_mc_result = 0.0;
		if(PRINT) print_mcdb_info(sim_params);
		if(SAVE_DATA) save_mcdb_data_fixed_tau(sim_params);
		sim_params.clear_mcdb_params();
		return;
	}


	while(sim_params.tau<MAX_TAU){

		sim_params.total_steps = MIN_STEPS_MCDB;
		sim_params.time_step = sim_params.tau/sim_params.total_steps;


		pre_exponentiate(sim_params);


		while(sim_params.total_steps <= MAX_STEPS_MCDB){

           // pre_exponentiate(sim_params);
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

			if(sim_params.new_distance > 0.3 and sim_params.total_steps >= MIN_STEPS_MCDB*4 and sim_params.total_steps <= MAX_STEPS_MCDB/2) break;
			sim_params.total_steps = sim_params.total_steps*2;
			sim_params.time_step = sim_params.tau/sim_params.total_steps;
			if(sim_params.total_steps <= MAX_STEPS_MCDB) scale_best_arrays_mcdb(sim_params, sim_params.j_best_scaled,sim_params.k_best_scaled,sim_params.b_best_scaled);

		}

		if(!update_distances(sim_params)) continue;
		if(PRINT) print_mc_results(sim_params);
		if(SAVE_DATA){
			sim_params.duration = (std::clock() - sim_params.start)/(double) CLOCKS_PER_SEC;
			save_mcdb_data_fixed_tau(sim_params);
		}
		if(sim_params.old_distance < sim_params.new_distance) break;
		if((sim_params.new_distance < DISTANCE_LIMIT && USE_ENERGY_DISTANCE) || (1-sim_params.evolved_target_dot_squared < DISTANCE_LIMIT && !USE_ENERGY_DISTANCE)) break;
		calc_tau(sim_params);
	}
	sim_params.clear_mcdb_params();
}




void mcdb_simulation(Simulation_Parameters& sim_params){
	int i,j, k,random_time_index, proposal_accepted,proposal_count;
	double *j_array, *k_array,*b_array,*j_temp, *k_temp, *b_temp, old_mc_result, new_mc_result, change;


	j_array          = new double[sim_params.total_steps]();
	k_array          = new double[sim_params.total_steps]();
	b_array          = new double[sim_params.total_steps]();
	j_temp           = new double[sim_params.total_steps]();
	k_temp           = new double[sim_params.total_steps]();
	b_temp           = new double[sim_params.total_steps]();
	sim_params.state = new double[2*sim_params.N]();



	if(sim_params.total_steps == MIN_STEPS_MCDB) init_arrays_mcdb(sim_params, sim_params.j_best,sim_params.k_best,sim_params.b_best);
	else copy_arrays_mcdb(sim_params,sim_params.j_best,sim_params.k_best,sim_params.b_best, sim_params.j_best_scaled,sim_params.k_best_scaled,sim_params.b_best_scaled,0,0);
	copy_arrays_mcdb(sim_params, j_array, k_array, b_array,sim_params.j_best,sim_params.k_best,sim_params.b_best,0,0);//a temporary array, used in the undoing of the changes

	std::memcpy(sim_params.state,sim_params.init_state, 2*sim_params.N*sizeof(double));
	evolve_mcdb(sim_params, j_array,k_array,b_array,0);
	sim_params.best_mc_result = cost(sim_params.target_state, sim_params.N, sim_params.state, sim_params.ham_target);
	old_mc_result = sim_params.best_mc_result;



	if(PRINT) print_mcdb_info(sim_params);


	for (i=0;i<TEMP_DECAY_ITERATIONS+ZERO_TEMP_ITERATIONS;i++){

		calc_new_temperature(sim_params, i);
		proposal_accepted = 0, proposal_count = 0;

		for (j=0; j<sim_params.total_sweeps*sim_params.sweeps_multiplier;j++){

			copy_arrays_mcdb(sim_params, j_temp, k_temp,b_temp,j_array, k_array, b_array,0,0);//a temporary array, used in the undoing of the changes

			random_time_index = change_array_mcdb(sim_params, j_array, k_array, b_array, j);

            std::memcpy(sim_params.state,sim_params.init_state, 2*sim_params.N*sizeof(double));
			evolve_mcdb(sim_params, j_array, k_array, b_array, 0);

			new_mc_result = cost(sim_params.target_state, sim_params.N,sim_params.state, sim_params.ham_target);

			if (new_mc_result<=sim_params.best_mc_result){
                sim_params.best_mc_result=new_mc_result;
                old_mc_result=new_mc_result;
                copy_arrays_mcdb(sim_params, sim_params.j_best, sim_params.k_best, sim_params.b_best,j_array, k_array, b_array, 0, 0);
                std::memcpy(sim_params.best_evolved_state,sim_params.state, 2*sim_params.N*sizeof(double));
            }
			else if (new_mc_result<=old_mc_result) old_mc_result=new_mc_result;
			else if (get_random_double(0,1,sim_params.rng)<exp(-(new_mc_result-old_mc_result)/(sim_params.temperature))) old_mc_result=new_mc_result, proposal_accepted++, proposal_count++;
			else copy_arrays_mcdb(sim_params, j_array, k_array, b_array,  j_temp, k_temp, b_temp, 0,0), proposal_count++;//undoing the change
		}
		if(PRINT) printf("           Best Cost:   %3.6f  ||  Acceptance Rate: %3.4f (%i/%i)\n",sim_params.best_mc_result,(double)proposal_accepted*1.0/proposal_count,proposal_accepted, proposal_count);
	}
	delete[] k_array, delete[] j_array, delete[] b_array, delete[] k_temp, delete[] j_temp, delete[] b_temp, delete[] sim_params.state;
}




void evolve_mcdb(Simulation_Parameters& sim_params, double *j_array, double *k_array, double *b_array, int random_time_index){
	int i=random_time_index, j, k,i_temp, x=0, s, index;
	double a, b, *pointer;
        
	while(i<sim_params.total_steps){
		i_temp = i;
		a = j_array[i], b = k_array[i];
		while(x+i<sim_params.total_steps){
			if(j_array[i+x] == a && k_array[i+x] == b) x+=1;
			else break;
		}
		i+=x;
		if(a != 0.0 || b != 0.0){
			if(a == 1.0 && b == 1.0) pointer = sim_params.e11;
			else if(a == 0.0 && b == 1.0) pointer = sim_params.e01;
			else if(a == 1.0 && b == 0.0) pointer = sim_params.e10;
			s = 1;
			while(x>0){
				while(s<=x) s = s*2;
				s = s/2;
				index = (int)round(log2( (double)(s*1.0*MAX_STEPS_MCDB/sim_params.total_steps))) * sim_params.N*sim_params.N*2;
				matrix_vector_mult(&pointer[index],sim_params.state, sim_params.N);
				x = x-s;
				s=1;
			}
		}
        x=0;
	}	
	if(CHECK) check_norm(sim_params.state, sim_params.N);
}




void calc_initial_temp_mcdb(Simulation_Parameters& sim_params){
	int i, j, random_time_index, count=0, k;
	double *state_random, *j_temp, *k_temp, *b_temp, *j_array,*k_array,*b_array,sum=0, change, old_mc_result, new_mc_result;

	j_array           = new double[sim_params.total_steps]();
	k_array           = new double[sim_params.total_steps]();
	b_array           = new double[sim_params.total_steps]();
	j_temp           = new double[sim_params.total_steps]();
	k_temp           = new double[sim_params.total_steps]();
	b_temp           = new double[sim_params.total_steps]();
	sim_params.state = new double[2*sim_params.N]();
	state_random     = new double[2*sim_params.N]();

	if(PRINT) printf("\n\n\n\n...Calculating initial temperature based on %i random starting states...\n", RANDOM_STATES);

	for (j=0;j<RANDOM_STATES;j++){
		for (i=0; i<sim_params.N*2;i++) state_random[i] = get_random_double(0,1,sim_params.rng);
        normalize_state(state_random, sim_params.N);
		std::memcpy(sim_params.state,state_random, 2*sim_params.N*sizeof(double));
		init_arrays_mcdb(sim_params, j_array, k_array, b_array);
		evolve_mcdb(sim_params,j_temp, k_temp, b_temp, 0);
		sim_params.best_mc_result = cost(sim_params.target_state, sim_params.N, sim_params.state, sim_params.ham_target);
        old_mc_result = sim_params.best_mc_result;

		for (i=0; i<sim_params.total_sweeps*sim_params.sweeps_multiplier;i++){

			copy_arrays_mcdb(sim_params, j_temp, k_temp,b_temp,j_array, k_array, b_array,0,0);//a temporary array, used in the undoing of the changes
		    random_time_index = change_array_mcdb(sim_params, j_array, k_array, b_array, i);

			std::memcpy(sim_params.state,state_random, 2*sim_params.N*sizeof(double));//resetting state
			evolve_mcdb(sim_params,j_array, k_array, b_array, 0);
			new_mc_result = cost(sim_params.target_state, sim_params.N, sim_params.state, sim_params.ham_target);


			if (new_mc_result<=sim_params.best_mc_result) sim_params.best_mc_result=new_mc_result, old_mc_result=new_mc_result;
			else if (new_mc_result<old_mc_result) old_mc_result=new_mc_result;
            else{
                sum += (new_mc_result-old_mc_result), count++;
                copy_arrays_mcdb(sim_params, j_array, k_array, b_array,  j_temp, k_temp, b_temp, 0,0);//undoing the change
            }
		}
	}
	sim_params.initial_temperature = -(sum/(count*log(ACCEPTANCE_PROB)));
	delete[] j_temp, delete[] k_temp, delete[] b_temp, delete[] sim_params.state, delete[] state_random;
}

\




void init_arrays_mcdb(Simulation_Parameters& sim_params, double *j_array,double *k_array,double *b_array){
	int i;
	for (i=0; i<sim_params.total_steps;i++){
		      j_array[i] = fmod(i+2.0, 2);//round(get_random_double(0,1.0,sim_params.rng));
	       	k_array[i] = fmod(i+1.0, 2);//round(get_random_double(0,1.0,sim_params.rng));
	       	b_array[i] = 1.0;//round(get_random_double(0,1.0,sim_params.rng));
	}
}




void copy_arrays_mcdb(Simulation_Parameters& sim_params, double* j_to,  double* k_to, double* b_to, double *j_from, double *k_from, double* b_from, int destination_index, int source_index){


	std::memcpy(&j_to[destination_index], &j_from[source_index], sizeof(double)*sim_params.total_steps);
	std::memcpy(&k_to[destination_index], &k_from[source_index], sizeof(double)*sim_params.total_steps);
	std::memcpy(&b_to[destination_index], &b_from[source_index], sizeof(double)*sim_params.total_steps);
}


void scale_best_arrays_mcdb(Simulation_Parameters& sim_params, double* j_best,double* k_best, double* b_best){
	int i;

	for(i=0;i<sim_params.total_steps/2;i++){
		j_best[(sim_params.total_steps-1)-2*i]     = j_best[(sim_params.total_steps/2-1)-i];
		j_best[(sim_params.total_steps-1)-(2*i+1)] = j_best[(sim_params.total_steps/2-1)-i];
		k_best[(sim_params.total_steps-1)-2*i]     = k_best[(sim_params.total_steps/2-1)-i];
		k_best[(sim_params.total_steps-1)-(2*i+1)] = k_best[(sim_params.total_steps/2-1)-i];
		b_best[(sim_params.total_steps-1)-2*i]     = b_best[(sim_params.total_steps/2-1)-i];
		b_best[(sim_params.total_steps-1)-(2*i+1)] = b_best[(sim_params.total_steps/2-1)-i];
	}
}


double change_array_mcdb(Simulation_Parameters& sim_params, double *j_array, double *k_array, double *b_array, int i){
	double *pointer, *pointer2;
	int random_time_index;

	random_time_index = (int)floor(get_random_double(0, sim_params.total_steps, sim_params.rng));


	if(i%2 == 0) pointer = j_array, pointer2 = k_array;
	if(i%2 == 1) pointer = k_array, pointer2 = j_array;
	//if(i%3 == 2) pointer = b_array;

	//rerolling the index if our current index has similar neighbors

	if(sim_params.total_steps == MAX_STEPS_MCDB and random_time_index > 0 and random_time_index < sim_params.total_steps-1){
		if(fmod(pointer[random_time_index] + pointer[random_time_index+1] + pointer[random_time_index-1], 3.0) <0.01){
			random_time_index = (int)floor(get_random_double(0, sim_params.total_steps, sim_params.rng));
		}
	}


	pointer[random_time_index] =  fmod(pointer[random_time_index] + 1.0, 2.0);
	//if(pointer[random_time_index]<0.01 and pointer2[random_time_index] <0.01) pointer[random_time_index] = 1;

	return random_time_index;
}



void pre_exponentiate(Simulation_Parameters& sim_params){
	int j;
	double *ham11,*ham01, *ham10,*ham_t_i,*jkb_vals;

	ham11  = new double[2*sim_params.N*sim_params.N]();
	ham10  = new double[2*sim_params.N*sim_params.N]();
	ham01  = new double[2*sim_params.N*sim_params.N]();
	ham_t_i      = new double[2*sim_params.N*sim_params.N]();
	jkb_vals     = new double[3]();


	for (j=0; j<sim_params.N*sim_params.N*2*TOTAL_STEP_CHANGES; j++) sim_params.e11[j]=0, sim_params.e01[j]=0,sim_params.e10[j]=0;




	int steps = MAX_STEPS_MCDB, index = 0;
	double time = 0;

	while(steps >= 1){
        jkb_vals[0] = 1, jkb_vals[1] = 1;
        construct_device_hamiltonian_uniform(sim_params, ham11, jkb_vals);
        jkb_vals[0] = 0, jkb_vals[1] = 1;
        construct_device_hamiltonian_uniform(sim_params, ham01, jkb_vals);
        jkb_vals[0] = 1, jkb_vals[1] = 0;
        construct_device_hamiltonian_uniform(sim_params, ham10, jkb_vals);
		
		time = sim_params.tau/steps;
		
    
		for (j=0; j<sim_params.N*sim_params.N; j++) ham_t_i[2*j+1] = (ham11[2*j]*-time), ham_t_i[2*j] = ham11[2*j+1]*time*-1; 
		exp_complex_double_matrix_pade(sim_params.N, ham_t_i, &sim_params.e11[index]);

		for (j=0; j<sim_params.N*sim_params.N; j++) ham_t_i[2*j+1] = (ham01[2*j]*-time), ham_t_i[2*j] = ham01[2*j+1]*time*-1;
		exp_complex_double_matrix_pade(sim_params.N, ham_t_i, &sim_params.e01[index]);

		for (j=0; j<sim_params.N*sim_params.N; j++) ham_t_i[2*j+1] = (ham10[2*j]*-time), ham_t_i[2*j] = ham10[2*j+1]*time*-1; 
		exp_complex_double_matrix_pade(sim_params.N, ham_t_i, &sim_params.e10[index]);

		if(CHECK){
			check_unitary(&sim_params.e11[index], sim_params.N);
			check_unitary(&sim_params.e01[index], sim_params.N);
			check_unitary(&sim_params.e10[index], sim_params.N);
		}

		steps = steps/2;
		index += 2*sim_params.N*sim_params.N;
	}
	delete[] jkb_vals, delete[] ham_t_i, delete[] ham11, delete[] ham10, delete[] ham01;
}


void get_sweeps_mcdb(Simulation_Parameters& sim_params){
	int x = sim_params.total_steps;
	double scalar, a = (double) MAX_STEPS_MCDB, b = (double) MIN_STEPS_MCDB, c = STEPS_CRUNCH_MCDB, offset;

	scalar = (a-b/c)/(a-b);
	offset = b/c - scalar*b;

	sim_params.total_sweeps = (int) ceil(scalar*x + offset) *SWEEPS_MCDB ;

}

