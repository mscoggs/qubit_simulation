#ifndef __MCBB_H_INCLUDED__
#define __MCBB_H_INCLUDED__

#include <gsl/gsl_rng.h>

#include "parameters.h"

/**
		Runs the monte carlo bang bang method, iterating over each seed. For each seed, the total time (tau) increased iteratively
        until the target ground state is met in the monte carlo bang bang simulation.

    @param sim_params contains all of the variables for the simulation
*/
void mcbb_method(Simulation_Parameters& sim_params);


/**
		Runs the monte carlo bang bang simluation for a fixed total time (tau). The temperature is annealed over each  iteration. This
        temperature dictactes the chance the a change in the j-k-b protocls (which drives us away from our desired state) is allowed.
        This process allows the protocls to get out of local minimum and hopefully find the global min.

    @param sim_params contains all of the variables for the simulation
*/
void mcbb_simulation(Simulation_Parameters& sim_params);


/**
    Evolves the starting state (the ground state of our initial hamiltonian) according to exp^(iH) where H is the hamiltonian which is
        a function of the j/k/b protocols

    @param sim_params contains all of the variables for the simulation
    @param j_array the protocol of the j parameter, where each odd-index holds the time of a jump (j=1) and the even-index holds the time of a drop(j=0)
    @param k_array the protocol of the k parameter, where each odd-index holds the time of a jump (k=1) and the even-index holds the time of a drop(k=0)
    @param b_array the protocol of the b parameter, where each odd-index holds the time of a jump (b=1) and the even-index holds the time of a drop(b=0)
*/
void evolve_mcbb(Simulation_Parameters& sim_params, double *j_array, double *k_array, double *b_array);


/**
    Calculates the initial temperature that will start the monte carlo simulation with roughly an acceptance rate of ACCEPTANCE_PROB (in parameters.h).

    @param sim_params contains all of the variables for the simulation
*/
void calc_initial_temp_mcbb(Simulation_Parameters& sim_params);


/**
    Once the ground state has been reached (within the limit of DIFFERENCE_LIMIT_MCBB), this function is called to search for the optimal time
        that achieves this.

    @param sim_params contains all of the variables for the simulation
*/
void binary_search_mcbb(Simulation_Parameters& sim_params);




/**
    Finds the next time where a change in one of the 3 protocols occurs given the current time

    @param time the current time in the evolution
    @param tau the total time allowed for the evolution
    @param j_array the protocol of the j parameter, where each odd-index holds the time of a jump (j=1) and the even-index holds the time of a drop(j=0)
    @param k_array the protocol of the k parameter, where each odd-index holds the time of a jump (k=1) and the even-index holds the time of a drop(k=0)
    @param b_array the protocol of the b parameter, where each odd-index holds the time of a jump (b=1) and the even-index holds the time of a drop(b=0)
    @param jkb_index the current index for each protocol (ie, the last time that was grabbed for each)
    @param jkb the values of the j/k/b arrays at a given time (either 1 or 0 for each, where jkb[0] is the j, jkb [1] is the k...)

    @return newtime the next time from the j/k/b arrays
*/
double find_next_time(double time, double tau,double* j_array, double* k_array,double*  b_array,int*  jkb_index, double*jkb);


/**
    Initializes the protocols to random values that are determined by the seed.

    @param sim_params contains all of the variables for the simulation
    @param j_array the protocol of the j parameter, where each odd-index holds the time of a jump (j=1) and the even-index holds the time of a drop(j=0)
    @param k_array the protocol of the k parameter, where each odd-index holds the time of a jump (k=1) and the even-index holds the time of a drop(k=0)
    @param b_array the protocol of the b parameter, where each odd-index holds the time of a jump (b=1) and the even-index holds the time of a drop(b=0)
*/
void init_arrays_mcbb(Simulation_Parameters& sim_params, double *j_array,double *k_array,double *b_array);


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
void copy_arrays_mcbb(Simulation_Parameters& sim_params, double* j_to,  double* k_to, double* b_to, double *j_from, double *k_from, double* b_from, int destination_index, int source_index);


/**
    Changes the arrays according randomly based on the random_time_index

    @param j_array the protocol of the j parameter, where each odd-index holds the time of a jump (j=1) and the even-index holds the time of a drop(j=0)
    @param k_array the protocol of the k parameter, where each odd-index holds the time of a jump (k=1) and the even-index holds the time of a drop(k=0)
    @param b_array the protocol of the b parameter, where each odd-index holds the time of a jump (b=1) and the even-index holds the time of a drop(b=0)
    @param change the amount the protocols will be changed by
    @param random_time_index the rng generated index which will pick the index of the array that gets changed (changing a row corresponds to changing all sites for a given timestep)
    @param i the iterator variable which cycles the type of change
    @param tau the total time allowed during evolution
*/
void change_array_mcbb(double *j_array, double *k_array, double *b_array, double change, int random_time_index, int i, double tau);


/**
    Calculates the change in time of the mcbb protocols, used in change_array_mcbb

    @param sim_params contains all of the variables for the simulation

    @return change the amount the protcol will be changed by
*/
double get_change_mcbb(Simulation_Parameters& sim_params);


/**
    Scales all three arrays, used when the total time is increased and we want to use an old jkb protocl as a inital guess for this new time

    @param j_to the protocol of the j parameter which will get written to, starting at destination_index
    @param k_to the protocol of the k parameter which will get written to, starting at destination_index
    @param b_to the protocol of the b parameter which will get written to, starting at destination_index
    @param scalar the amount that each time will be scaled by
*/
void scale_arrays_mcbb(double *j_array, double *k_array, double *b_array, double scalar);


/**
	prints tau vs energy for the evolution of a fixed protocol, where the time of the jumps are a fixed fraction of tau
    @param sim_params contains all of the variables for the simulation
*/
void evolve_fixed_protocol(Simulation_Parameters& sim_params);

#endif
