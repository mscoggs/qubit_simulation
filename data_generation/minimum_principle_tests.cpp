#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <string>
#include <algorithm>
#include <ctime>

#include "check.h"
#include "hamiltonian.h"
#include "operations.h"
#include "parameters.h"
#include "print.h"
#include "linear_algebra.h"



main (int argc, char *argv[]){
    double tau = 0.351562, time_step = 0.00001, time = tau, re =0, im = 0;
    int i,j,j_index = 0, k_index=3, total_steps = (int)floor(tau/time_step), N = 36;
    int h_index = total_steps-1;
    double *H_j, *H_k, *pi, *psi, *psi_temp1, *psi_temp2, *hamiltonian,*ham_t_i, *ham_real,*exp_matrix, *hamiltonian10, *hamiltonian01;
    double target_state[N*2] = {-0.0392708, 0, -0.0354592, 0, -0.0392708, 0, -0.136812, 0, -0.244756, 0, -0.0392708, 0, -0.244756, 0, -0.198235, 0, -0.244756, 0, -0.198235, 0, -0.0392708, 0, -0.244756, 0, -0.136812, 0, -0.308535, 0, -0.198235, 0, -0.0354592, 0, -0.136812, 0, -0.0551019, 0, -0.136812, 0, -0.244756, 0, -0.0392708, 0, -0.136812, 0, -0.308535, 0, -0.136812, 0, -0.244756, 0, -0.198235, 0, -0.244756, 0, -0.0392708, 0, -0.0551019, 0, -0.136812, 0, -0.0354592, 0, -0.0392708, 0, -0.244756, 0, -0.136812, 0, -0.0354592, 0, -0.0392708, 0};
    double evolved_state[N*2] = {0.0050459, 0.0358555, 0.0492864, 0.0336511, 0.0050459, 0.0358555, 0.0505893, 0.153499, 0.0647756, 0.227186, 0.0050459, 0.0358555, 0.0647756, 0.227186, 0.0234769, 0.189666, 0.0647756, 0.227186, 0.0234769, 0.189666, 0.0050459, 0.0358555, 0.0647756, 0.227186, 0.0505893, 0.153499, 0.0533961, 0.279351, 0.0234769, 0.189666, 0.0492864, 0.0336511, 0.0505893, 0.153499, 0.0449853, 0.062841, 0.0505893, 0.153499, 0.0647756, 0.227186, 0.0050459, 0.0358555, 0.0505893, 0.153499, 0.0533961, 0.279351, 0.0505893, 0.153499, 0.0647756, 0.227186, 0.0234769, 0.189666, 0.0647756, 0.227186, 0.0050459, 0.0358555, 0.0449853, 0.062841, 0.0505893, 0.153499, 0.0492864, 0.0336511, 0.0050459, 0.0358555, 0.0647756, 0.227186, 0.0505893, 0.153499, 0.0492864, 0.0336511, 0.0050459, 0.0358555};
    double k_on[3]     = {0,1,0};
    double jkb_vals[3] = {1,0,0};
    double j_opt[1] = {0.00948312};
    double k_opt[4] = {0, 0.0983751, 0.134202, 0.316601};
    double j_on[3]     = {1,0,0};


    hamiltonian = new double[2*sim_params.N*sim_params.N]();
    hamiltonian10 = new double[2*sim_params.N*sim_params.N]();
    hamiltonian01 = new double[2*sim_params.N*sim_params.N]();
    exp_matrix  = new double[2*sim_params.N*sim_params.N]();
    ham_t_i     = new double[2*sim_params.N*sim_params.N]();
    ham_real    = new double[sim_params.N*sim_params.N]();
    H_j           = new double[total_steps]();
    H_k           = new double[total_steps]();
    pi           = new double[2*N]();
    psi           = new double[2*N]();
    psi_temp1           = new double[2*N]();
    psi_temp2           = new double[2*N]();

    construct_device_hamiltonian_uniform(sim_params, hamiltonian10, j_on);
    construct_device_hamiltonian_uniform(sim_params, hamiltonian01, k_on);

    Simulation_Parameters sim_params;
    sim_params.initialize_lattice(2,sim_params);

    for(i=0;i<N*2;i++) psi[i] = evolved_state[i];

    for(i=0;i<N*2;i+=2){
            re += target_state[i]*evolved_state[i] + target_state[i+1]*evolved_state[i+1];
            im += target_state[i]*evolved_state[i+1] - target_state[i+1]*evolved_state[i];
    }
    for(i=0;i<N*2;i+=2){
        pi[i] = -2 * (re*target_state[i] - im*target_state[i+1]);
        pi[i+1] = -2 * (re*target_state[i+1] + im*target_state[i]);
    }

    while(time >= 0){
        ///manage jkb_vals
        if (j_index >= 0) {
            if(time < j_opt[j_index]){
                jkb_vals[0] = fmod(jkb_vals[0] +1, 2.0);
                j_index -= 1;
            }
        }
        if (k_index >= 0) {
            if(time < k_opt[k_index]){
                jkb_vals[1] = fmod(jkb_vals[1] +1, 2.0);
                k_index -= 1;
            }
        }

        construct_device_hamiltonian_uniform(sim_params, hamiltonian, jkb_vals);
        for (j=0; j<sim_params.N*sim_params.N; j++) ham_t_i[2*j+1] = (hamiltonian[2*j]*time_step), ham_t_i[2*j] = -(hamiltonian[2*j+1]*time_step);  //multiplying by i*dt for the Pade approximation
        exp_complex_double_matrix_pade(sim_params.N, ham_t_i, exp_matrix);
        if(CHECK) check_unitary(exp_matrix, sim_params.N);

        matrix_vector_mult(exp_matrix, pi, sim_params.N);
        matrix_vector_mult(exp_matrix, psi, sim_params.N);


        memcpy(psi_temp1, psi, sizeof(double)*2*N);
        memcpy(psi_temp2, psi, sizeof(double)*2*N);

        if(CHECK) check_norm(psi, N);

        construct_device_hamiltonian_uniform(sim_params, hamiltonian10, j_on);
        construct_device_hamiltonian_uniform(sim_params, hamiltonian01, k_on);

        matrix_vector_mult(hamiltonian10, psi_temp1, sim_params.N);
        matrix_vector_mult(hamiltonian01, psi_temp2, sim_params.N);

        for (i=0; i<2*N; i+=2){
            H_j[h_index] += pi[i] * (psi_temp1[i+1])   - pi[i+1] * (psi_temp1[i]);
            H_k[h_index] += pi[i] * (psi_temp2[i+1])   - pi[i+1] * (psi_temp2[i]);
        }

        h_index -= 1;
        time  -= time_step;
    }


    std::ofstream file;
    file.open("../data_paper_98/h_data.txt");
    file << "H_J = [";
    for (j=0; j<total_steps; j++) file << H_j[j] << ", ";

    file << "]\nH_K = [";
    for (j=0; j<total_steps; j++) file << H_k[j] << ", ";
    file << "]";
}
