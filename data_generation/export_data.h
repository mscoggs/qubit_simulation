#ifndef __EXPORT_H_INCLUDED__
#define __EXPORT_H_INCLUDED__

#include <string>

#include "parameters.h"


/**
    Saves the monte carlo brute force simulation data

    @param sim_params contains all of the variables for the simulation
*/
void save_mcbf_data(Simulation_Parameters& sim_params);


/**
    Saves the monte carlo bang bang simulation data

    @param sim_params contains all of the variables for the simulation
*/
void save_mcbb_data(Simulation_Parameters& sim_params);


/**
    Saves the adiabatic simulation data

    @param sim_params contains all of the variables for the simulation
*/
void save_adiabatic_data(Simulation_Parameters& sim_params);


/**
    Generates a path to save the data that depends on the type of simulation and the jkb values of the start and target hamiltonians

    @param sim_params contains all of the variables for the simulation
    @param type the type of simulation
*/
std::string make_path(Simulation_Parameters sim_params, std::string type);


/**
    A sub function of all the save_blank_data functions, which saves the start and target jkb values along with some info about the lattice

    @param sim_params contains all of the variables for the simulation
    @param path the location of the file that will be appended with the ham parameters
*/
void save_hamiltonian_parameters(Simulation_Parameters sim_params,std::string path);


#endif
