#ifndef __MCDB_H_INCLUDED__
#define __MCDB_H_INCLUDED__

#include <gsl/gsl_rng.h>

#include "parameters.h"

/**
		Runs the monte carlo discrete bang method, iterating over each seed. For each seed, the total time (tau) increased iteratively
        until the target ground state is met in the monte carlo bang bang simulation.

    @param sim_params contains all of the variables for the simulation
*/
void mcdb_method(Simulation_Parameters& sim_params);


/**
		Runs the monte carlo bang bang simluation for a fixed total time (tau). The temperature is annealed over each  iteration. This
        temperature dictactes the chance the a change in the j-k-b protocls (which drives us away from our desired state) is allowed.
        This process allows the protocls to get out of local minimum and hopefully find the global min.

    @param sim_params contains all of the variables for the simulation
*/
void mcdb_simulation(Simulation_Parameters& sim_params);


/**
    Evolves the starting state (the ground state of our initial hamiltonian) according to exp^(iH) where H is the hamiltonian which is
        a function of the j/k/b protocols

    @param sim_params contains all of the variables for the simulation
    @param j_array the protocol of the j parameter, where each odd-index holds the time of a jump (j=1) and the even-index holds the time of a drop(j=0)
    @param k_array the protocol of the k parameter, where each odd-index holds the time of a jump (k=1) and the even-index holds the time of a drop(k=0)
    @param b_array the protocol of the b parameter, where each odd-index holds the time of a jump (b=1) and the even-index holds the time of a drop(b=0)
    @param random_time_index this is the index that was changed, and we skip the evolution prior to the change by pulling from the saved_states
*/
void evolve_mcdb(Simulation_Parameters& sim_params, double *j_array, double *k_array, double *b_array, int random_time_index);


/**
    Calculates the initial temperature that will start the monte carlo simulation with roughly an acceptance rate of ACCEPTANCE_PROB_MCBB (in parameters.h).

    @param sim_params contains all of the variables for the simulation
*/
void calc_initial_temp_mcdb(Simulation_Parameters& sim_params);



/**
    Calculates the tau for the next mcdb_simulation

    @param sim_params contains all of the variables for the simulation
*/
void calc_tau_mcdb(Simulation_Parameters& sim_params);




/**
    Initializes the protocols to random values that are determined by the seed.

    @param sim_params contains all of the variables for the simulation
    @param j_array the protocol of the j parameter, where each odd-index holds the time of a jump (j=1) and the even-index holds the time of a drop(j=0)
    @param k_array the protocol of the k parameter, where each odd-index holds the time of a jump (k=1) and the even-index holds the time of a drop(k=0)
    @param b_array the protocol of the b parameter, where each odd-index holds the time of a jump (b=1) and the even-index holds the time of a drop(b=0)
*/
void init_arrays_mcdb(Simulation_Parameters& sim_params, double *j_array,double *k_array,double *b_array);


/**
    Copies the info from a set of j/k/b arrays to another set, allowing for an offset to and from each

    @param sim_params contains all of the variables for the simulation
    @param j_to the protocol of the j parameter which will get written to, starting at destination_index
    @param k_to the protocol of the k parameter which will get written to, starting at destination_index
    @param b_to the protocol of the b parameter which will get written to, starting at destination_index
    @param j_from the protocol of the j parameter which will get written from, starting at 2source_index
    @param k_from the protocol of the k parameter which will get written from, starting at source_index
    @param b_from the protocol of the b parameter which will get written from, starting at source_index
    @param destination_index the offset for the arrays that are getting copied to
    @param source_index the offset for the arrays that are getting copied from
*/
void copy_arrays_mcdb(Simulation_Parameters& sim_params, double* j_to,  double* k_to, double* b_to, double *j_from, double *k_from, double* b_from, int destination_index, int source_index);


/**
    scale the arrays, cutting each block in half and maintaining the previous values

    @param sim_params contains all of the variables for the simulation
*/
void scale_best_arrays_mcdb(Simulation_Parameters& sim_params, double* j_best,  double* k_best, double* b_best);




/**
    Changes the arrays according randomly based on the random_time_index

    @param sim_params contains all of the variables for the simulation
    @param j_array the protocol of the j parameter, where each odd-index holds the time of a jump (j=1) and the even-index holds the time of a drop(j=0)
    @param k_array the protocol of the k parameter, where each odd-index holds the time of a jump (k=1) and the even-index holds the time of a drop(k=0)
    @param b_array the protocol of the b parameter, where each odd-index holds the time of a jump (b=1) and the even-index holds the time of a drop(b=0)
    @param i the iterator variable which cycles the type of change
*/
double change_array_mcdb(Simulation_Parameters& sim_params, double *j_array, double *k_array, double *b_array, int i);



/**
    loads up the exponentiate matrices

    @param sim_params contains all of the variables for the simulation

*/
void pre_exponentiate(Simulation_Parameters& sim_params);


/**
    scales the total sweeps, 'crunching' them in the beginning, then increasing them as total steps increase

    @param sim_params contains all of the variables for the simulation

*/
void get_sweeps_mcdb(Simulation_Parameters& sim_params);
#endif
