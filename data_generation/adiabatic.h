#ifndef __ADIA_H_INCLUDED__
#define __ADIA_H_INCLUDED__

#include <gsl/gsl_rng.h>

#include "parameters.h"


/**
    Runs the adiabatic method which iteratively increases total time (tau) until the ground state is reached after evolution, using a linear parameterization

    @param sim_params contains all of the variables for the simulation
*/
void adiabatic_method(Simulation_Parameters& sim_params);


/**
    Evolves the starting state (the ground state of our initial hamiltonian) according to exp^(iH) where H is the hamiltonian which is
        a function of the j/k/b protocols. In the adiabatic evolution, the values of j/k/b are slowly changed from the initial to the target values

    @param sim_params contains all of the variables for the simulation
*/
void evolve_adiabatic(Simulation_Parameters& sim_params);

#endif
