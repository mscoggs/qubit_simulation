#ifndef __PRINT_H_INCLUDED__
#define __PRINT_H_INCLUDED__

#include "parameters.h"


/**
   Prints an N dimensional state represented in an 2*N array, with odd indices holding the imaginary components and even holding the real components
	 Ex: state[2n] = 1, state[2n+1]=4 means the nth entry of this state is 1+4i.

   @param state the quantum state to be printed
   @param N the dimension of the state
 */
void print_state(double* state,int N);


/**
    Prints a complex hamiltonian in matrix form with the odd indices holding the complex components and the even holding the real components
	  Ex: hamiltonian[2*N*C + 2*R] = 2, hamiltonian[2*N*C + 2*R  + 1] = 4 means the entry at column=C and row=R is 2+4i, with 0 indexing on C and R

    @param hamiltonian the 2*Nx2*N complex hamiltonian to be printed
    @param N the dimension of the quantum system
 */
void print_hamiltonian_complex(double* hamiltonian, int N);


/**
    Prints a real hamiltonian in matrix form where hamiltonian[N*C + R] accesses the Cth column and the Rth row (using 0 indexing)

    @param hamiltonian the NxN real hamiltonian to be printed
    @param N the dimension of the quantum system
 */
void print_hamiltonian_real(double *hamiltonian, int N);


/**
    Prints the mcbb arrays, where each index holds a time that the value of j/k/b will change.

    @param j_times the array holding the time of each jump/drop for the j parameter. The even indicies hold a jump, the odds hold a drop.
    @param k_times the array holding the time of each jump/drop for the k parameter. The even indicies hold a jump, the odds hold a drop.
    @param b_times the array holding the time of each jump/drop for the b parameter. The even indicies hold a jump, the odds hold a drop.
 */
void print_arrays_mcbb(double* j_times, double* k_times, double* b_times);

/**
    Prints the mcdb arrays, where each index holds the value of j/k/b during evolution

    @param j_times the array holding the time of each jump/drop for the j parameter. The even indicies hold a jump, the odds hold a drop.
    @param k_times the array holding the time of each jump/drop for the k parameter. The even indicies hold a jump, the odds hold a drop.
    @param b_times the array holding the time of each jump/drop for the b parameter. The even indicies hold a jump, the odds hold a drop.
    @oaram size the size of the arrays
 */
void print_arrays_mcdb(double* j_times, double* k_times, double* b_times, int size);


/**
    Prints the mcbf arrays, where each index holds a value at a certain time of the evolution.

    @param j_times the array holding the value of each step for the j parameter.
    @param k_times the array holding the value of each step for the k parameter.
    @param b_times the array holding the value of each step for the b parameter.
    @param total_steps the total number of steps for the evolution, the total number of inputs for each array is equal to
           total_steps*NUMBER_OF_SITES for b (a site dependent parameter), and twice that for j/k, a bond dependent parameter
 */
void print_arrays_mcbf(double *j_array, double *k_array, double *b_array, int total_step);


/**
    Prints the results of the monte-carlo run

    @param sim_params contains all of the variables for the simulation
 */
void print_mc_results(Simulation_Parameters& sim_params);


/**
    Prints the info relevant to the adiabatic simulation

    @param sim_params contains all of the variables for the simulation
 */
void print_adiabatic_info(Simulation_Parameters& sim_params);


/**
    Prints the info relevant to the monte-carlo brute-force simulation

    @param sim_params contains all of the variables for the simulation
 */
void print_mcbf_info(Simulation_Parameters& sim_params);


/**
    Prints the info relevant to the monte-carlo discrete-bang simulation

    @param sim_params contains all of the variables for the simulation
 */
void print_mcdb_info(Simulation_Parameters& sim_params);


/**
    Prints the info relevant to the monte-carlo bang-bang simulation

    @param sim_params contains all of the variables for the simulation
 */
void print_mcbb_info(Simulation_Parameters& sim_params);

#endif
