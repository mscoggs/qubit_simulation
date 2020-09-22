#ifndef __MCBF_H_INCLUDED__
#define __MCBF_H_INCLUDED__

#include <gsl/gsl_rng.h>

#include "parameters.h"


/**
		Runs the monte carlo brute force method, iterating over each seed. For each seed, the total time (tau) increased iteratively
        until the target ground state is met in the monte carlo brute force simulation.

    @param sim_params contains all of the variables for the simulation
*/
void mcbf_method(Simulation_Parameters& sim_params);

/**
		Runs the monte carlo brute force simluation for a fixed total time (tau). The temperature is annealed over each  iteration. This
        temperature dictactes the chance the a change in the j-k-b protocls (which drives us away from our desired state) is allowed.
        This process allows the protocls to get out of local minimum and hopefully find the global min.

    @param sim_params contains all of the variables for the simulation
*/
void mcbf_simulation(Simulation_Parameters& sim_params);


/**
    Evolves the starting state (the ground state of our initial hamiltonian) according to exp^(iH) where H is the hamiltonian which is
        a function of the j/k/b protocols

    @param sim_params contains all of the variables for the simulation
    @param j_array the protocol of the j parameter, where each element holds the value of the parameter at a given time during the evolution
    @param k_array the protocol of the k parameter, where each element holds the value of the parameter at a given time during the evolution
    @param b_array the protocol of the b parameter, where each element holds the value of the parameter at a given time during the evolution
*/
void evolve_mcbf(Simulation_Parameters& sim_params, double *j_array, double *k_array, double *b_array);


/**
    Calculates the initial temperature that will start the monte carlo simulation with roughly an acceptance rate of ACCEPTANCE_PROB (in parameters.h).

    @param sim_params contains all of the variables for the simulation
*/
void calc_initial_temp_mcbf(Simulation_Parameters& sim_params);


/**
    Once the ground state has been reached (within the limit of DIFFERENCE_LIMIT_MC), this function is called to search for the optimal time
        that achieves this.

    @param sim_params contains all of the variables for the simulation
*/
void binary_search_mcbf(Simulation_Parameters& sim_params);


/**
    Calculates the total_steps for the next mcbf_simulation

    @param sim_params contains all of the variables for the simulation


*/
void calc_total_steps_mcbf(Simulation_Parameters& sim_params);


/**
    Initializes the protocols to random values that are determined by the seed.

    @param sim_params contains all of the variables for the simulation
    @param j_array the protocol of the j parameter, where each element holds the value of the parameter at a given time during the evolution
    @param k_array the protocol of the k parameter, where each element holds the value of the parameter at a given time during the evolution
    @param b_array the protocol of the b parameter, where each element holds the value of the parameter at a given time during the evolution
    @param total_steps the total number of steps for the evolution, the total number of inputs for each array is equal to
           total_steps*NUMBER_OF_SITES for b (a site dependent parameter), and twice that for j/k, a bond dependent parameter
*/
void init_arrays_mcbf(Simulation_Parameters& sim_params, double *j_array, double *k_array, double *b_array);


/**
    Copies the info from a set of j/k/b arrays to another set, allowing for an offset to and from each

    @param sim_params contains all of the variables for the simulation
    @param j_to the protocol of the j parameter which will get written to, starting at 2*destination_index (a factor of 2 because the j/k arrays are twice as big)
    @param k_to the protocol of the k parameter which will get written to, starting at 2*destination_index
    @param b_to the protocol of the b parameter which will get written to, starting at destination_index
    @param j_from the protocol of the j parameter which will get written from, starting at 2*source_index
    @param k_from the protocol of the k parameter which will get written from, starting at 2*source_index
    @param b_from the protocol of the b parameter which will get written from, starting at source_index
    @param destination_index the offset for the arrays that are getting copied to
    @param source_index the offset for the arrays that are getting copied from
*/
void copy_arrays_mcbf(Simulation_Parameters& sim_params, double* j_to,  double* k_to, double* b_to, double *j_from, double *k_from, double* b_from, int destination_index, int source_index);


/**
    Changes the arrays according to the constants set in parameters.h

    @param j_array the protocol of the j parameter, where each element holds the value of the parameter at a given time during the evolution
    @param k_array the protocol of the k parameter, where each element holds the value of the parameter at a given time during the evolution
    @param b_array the protocol of the b parameter, where each element holds the value of the parameter at a given time during the evolution
    @param random_row the rng generated row which will get changed (changing a row corresponds to changing all sites for a given timestep)
    @param random_col the rng generated column which will get changed (changing a column corresponds to changing all timesteps for a given site)
    @param change the amount the protocols will be changed by
    @param i the iterator variable which cycles the type of change
    @param total_steps the total number of steps for the evolution, the total number of inputs for each array is equal to
           total_steps*NUMBER_OF_SITES for b (a site dependent parameter), and twice that for j/k, a bond dependent parameter
*/
void change_array_mcbf(double *j_array, double *k_array, double *b_array, int random_row, int random_col, double change, int i, int total_steps);


/**
    Changes the rows of the protocols, where a row corresponds to all sites for a given timestep. This function is called by change_array_mcbf
        according to the parameters set in parameters.h under the mcbf section

    @param j_array the protocol of the j parameter, where each element holds the value of the parameter at a given time during the evolution
    @param k_array the protocol of the k parameter, where each element holds the value of the parameter at a given time during the evolution
    @param b_array the protocol of the b parameter, where each element holds the value of the parameter at a given time during the evolution
    @param row the row (timestep) that will be changed
    @param change the amount the protocols will be changed by
    @param k a bool that determines if the k array will be changed
    @param k a bool that determines if the k array will be changed
    @param k a bool that determines if the k array will be changed
    @param jump the number of elements that will be jumped over, if set to 1 this function changes all in a given row. If 3, it changes every 3rd, etc.
    @param offset the starting location of the change. Ex: If this is 0 and jump is 2, all even elements in row will be changed,
        if this is 1 and jump is 2, all odd elements will be changed
*/
void change_row_mcbf(double *j_array, double *k_array,double *b_array, int row, double change, bool j, bool k, bool b, int jump, int offset);


/**
    Changes the columns of the protocols, where a column corresponds to all timesteps for a given site. This function is called by change_array_mcbf
        according to the parameters set in parameters.h under the mcbf section

    @param j_array the protocol of the j parameter, where each element holds the value of the parameter at a given time during the evolution
    @param k_array the protocol of the k parameter, where each element holds the value of the parameter at a given time during the evolution
    @param b_array the protocol of the b parameter, where each element holds the value of the parameter at a given time during the evolution
    @param col the col (site) that will be changed
    @param change the amount the protocols will be changed by
    @param k a bool that determines if the k array will be changed
    @param k a bool that determines if the k array will be changed
    @param k a bool that determines if the k array will be changed
    @param jump the number of elements that will be jumped over, if set to 1 this function changes all in a given column. If 3, it changes every 3rd, etc.
    @param offset the starting location of the change. Ex: If this is 0 and jump is 2, all even elements in column will be changed,
        if this is 1 and jump is 2, all odd elements will be changed
*/
void change_col_mcbf(int total_steps,double *j_array,double *k_array,double *b_array, int col, double change,bool j, bool k, bool b, int jump, int offset);


/**
    Changes all 3 arrays for a single time and a single site

    @param j_array the protocol of the j parameter, where each element holds the value of the parameter at a given time during the evolution
    @param k_array the protocol of the k parameter, where each element holds the value of the parameter at a given time during the evolution
    @param b_array the protocol of the b parameter, where each element holds the value of the parameter at a given time during the evolution
    @param row the row (timestep) that will be changed
    @param col the col (site) that will be changed
    @param change the amount the protocols will be changed by
*/
void change_single_mcbf(double *j_array,double *k_array,double *b_array, int row,int col, double change);


/**
    Calculates the change in the mcbf protocols, used in change_array_mcbf

    @param sim_params contains all of the variables for the simulation

    @return change the amount the protcol will be changed by
*/
double get_change_mcbf(Simulation_Parameters& sim_params, int temp_index);


/**
    Scales the protocols, allowing for a finer discretization

    @param sim_params contains all of the variables for the simulation
    @param j_array the protocol of the j parameter, where each element holds the value of the parameter at a given time during the evolution
    @param k_array the protocol of the k parameter, where each element holds the value of the parameter at a given time during the evolution
    @param b_array the protocol of the b parameter, where each element holds the value of the parameter at a given time during the evolution
    @param total_steps the total number of steps for the evolution, the total number of inputs for each array is equal to
           total_steps*NUMBER_OF_SITES for b (a site dependent parameter), and twice that for j/k, a bond dependent parameter
*/
void scale_arrays_mcbf(double* j_array, double* k_array,double* b_array, int total_steps);


#endif
