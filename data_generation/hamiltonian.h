#ifndef __HAM_H_INCLUDED__
#define __HAM_H_INCLUDED__


#include "parameters.h"


/**
    Generating the XXZ hamiltonian allowing for non-uniform sites for j/k/b
    H_dev = sum_ij {  j * (S_i^x * S_j^x  +  S_i^y * S_j^y)  + k * (S_i^z * S_j^z) } + sum_i{b * S_i^Z}

    @param sim_params contains all of the variables for the simulation
    @param hamiltonian the array that will store the hamiltonian
    @param j_array the protocol of the j parameter, where each element holds the value of the parameter at a given time during the evolution
    @param k_array the protocol of the k parameter, where each element holds the value of the parameter at a given time during the evolution
    @param b_array the protocol of the b parameter, where each element holds the value of the parameter at a given time during the evolution
    @param index the index that we'll be using in the j/k/b arrays for the generation of the hamiltonian.
*/
void construct_device_hamiltonian(Simulation_Parameters sim_params, double *hamiltonian ,double *j_array, double *k_array, double *b_array,  int index);


/**
    Generating the XXZ hamiltonian for uniform sites for j/k/b
    H_dev = sum_ij {  j * (S_i^x * S_j^x  +  S_i^y * S_j^y)  + k * (S_i^z * S_j^z) } + sum_i{b * S_i^Z}

    @param sim_params contains all of the variables for the simulation
    @param hamiltonian the array that will store the hamiltonian
    @param jkb the values of j, k, and b for the hamiltonian where (jkb[0] = j, jkb[1] = k, ...)
*/
void construct_device_hamiltonian_uniform(Simulation_Parameters sim_params, double *hamiltonian , double *jkb);


/**
    Generating the Model hamiltonian for uniform sites for parameters T and V
    H_model = sum_ij {  -t * (c_i^dagger * c_j  + c_j^dagger * c_i )  + v * (n_i - .5)(n_j - .5 }
    where a Jordan-Wigner transformation is used to map these to spins

    @param sim_params contains all of the variables for the simulation
    @param hamiltonian the array that will store the hamiltonian
*/
void construct_model_hamiltonian(Simulation_Parameters sim_params, double *hamiltonian);

#endif
