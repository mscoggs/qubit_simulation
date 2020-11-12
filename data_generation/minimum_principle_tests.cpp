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
  double tau = 0.351562, time_step = 0.00001, time = tau;
  int i,j;
	double re =0, im = 0;
  int total_steps = (int)floor(tau/time_step), N = 36;
  double *H_j, *H_k, *pi, *psi, *psi_temp1, *psi_temp2;
  double target_state[N*2] = {-0.0392708, 0, -0.0354592, 0, -0.0392708, 0, -0.136812, 0, -0.244756, 0, -0.0392708, 0, -0.244756, 0, -0.198235, 0, -0.244756, 0, -0.198235, 0, -0.0392708, 0, -0.244756, 0, -0.136812, 0, -0.308535, 0, -0.198235, 0, -0.0354592, 0, -0.136812, 0, -0.0551019, 0, -0.136812, 0, -0.244756, 0, -0.0392708, 0, -0.136812, 0, -0.308535, 0, -0.136812, 0, -0.244756, 0, -0.198235, 0, -0.244756, 0, -0.0392708, 0, -0.0551019, 0, -0.136812, 0, -0.0354592, 0, -0.0392708, 0, -0.244756, 0, -0.136812, 0, -0.0354592, 0, -0.0392708, 0};
  double evolved_state[N*2] = {0.0050459, 0.0358555, 0.0492864, 0.0336511, 0.0050459, 0.0358555, 0.0505893, 0.153499, 0.0647756, 0.227186, 0.0050459, 0.0358555, 0.0647756, 0.227186, 0.0234769, 0.189666, 0.0647756, 0.227186, 0.0234769, 0.189666, 0.0050459, 0.0358555, 0.0647756, 0.227186, 0.0505893, 0.153499, 0.0533961, 0.279351, 0.0234769, 0.189666, 0.0492864, 0.0336511, 0.0505893, 0.153499, 0.0449853, 0.062841, 0.0505893, 0.153499, 0.0647756, 0.227186, 0.0050459, 0.0358555, 0.0505893, 0.153499, 0.0533961, 0.279351, 0.0505893, 0.153499, 0.0647756, 0.227186, 0.0234769, 0.189666, 0.0647756, 0.227186, 0.0050459, 0.0358555, 0.0449853, 0.062841, 0.0505893, 0.153499, 0.0492864, 0.0336511, 0.0050459, 0.0358555, 0.0647756, 0.227186, 0.0505893, 0.153499, 0.0492864, 0.0336511, 0.0050459, 0.0358555};
  H_j           = new double[total_steps]();
  H_k           = new double[total_steps]();
  pi           = new double[2*N]();
  psi           = new double[2*N]();
  psi_temp1           = new double[2*N]();
  psi_temp2           = new double[2*N]();



  for(i=0;i<N*2;i++){
    psi[i] = evolved_state[i];
  }

	for(i=0;i<N*2;i+=2){
		re += target_state[i]*evolved_state[i] + target_state[i+1]*evolved_state[i+1];
		im += target_state[i]*evolved_state[i+1] - target_state[i+1]*evolved_state[i];
	}
  for(i=0;i<N*2;i+=2){
    pi[i] = -2 * (re*target_state[i] - im*target_state[i+1]);
    pi[i+1] = -2 * (re*target_state[i+1] + im*target_state[i]);

  }

  Simulation_Parameters sim_params;
	sim_params.initialize_lattice(2,sim_params);





  double *hamiltonian,*ham_t_i, *ham_real,*exp_matrix, *hamiltonian10, *hamiltonian01;

  hamiltonian = new double[2*sim_params.N*sim_params.N]();
  hamiltonian10 = new double[2*sim_params.N*sim_params.N]();
  hamiltonian01 = new double[2*sim_params.N*sim_params.N]();
	exp_matrix  = new double[2*sim_params.N*sim_params.N]();
	ham_t_i     = new double[2*sim_params.N*sim_params.N]();
	ham_real    = new double[sim_params.N*sim_params.N]();


  double jkb_vals[3] = {1,0,0};
  double j_opt[1] = {0.00948312};
  double k_opt[4] = {0, 0.0983751, 0.134202, 0.316601};
  double j_on[3]     = {1,0,0};
  double k_on[3]     = {0,1,0};
  int j_index = 0;
  int k_index = 3;

  int h_index = total_steps-1;

  construct_device_hamiltonian_uniform(sim_params, hamiltonian10, j_on);
  construct_device_hamiltonian_uniform(sim_params, hamiltonian01, k_on);


  while(time >= 0){
    ///manage jkb_vals
    if (j_index >= 0) {
      if(time < j_opt[j_index]){
        jkb_vals[0] = fmod(jkb_vals[0] +1, 2.0);
        j_index -= 1;
        //printf("time : %f, j: %f, k: %f\n", time, jkb_vals[0], jkb_vals[1]);
      }
    }
    if (k_index >= 0) {
      if(time < k_opt[k_index]){
        jkb_vals[1] = fmod(jkb_vals[1] +1, 2.0);
        k_index -= 1;
        //printf("time : %f, j: %f, k: %f\n", time, jkb_vals[0], jkb_vals[1]);
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
    // print_state(pi, N);
    // print_state(psi, N);
    //
    // exit(0);


    for (i=0; i<2*N; i+=2){
      // double j_on_print = pi[i] * (psi_temp1[i+1])   - pi[i+1] * (psi_temp1[i]);
      // printf("jon = %f\n", j_on_print);
      // //if(i == 8) exit(0);
      H_j[h_index] += pi[i] * (psi_temp1[i+1])   - pi[i+1] * (psi_temp1[i]);
      H_k[h_index] += pi[i] * (psi_temp2[i+1])   - pi[i+1] * (psi_temp2[i]);
    }

    h_index -= 1;
    time  -= time_step;
  }
  print_state(psi, N);




  std::ofstream file;
	file.open("../data_paper_98/h_data.txt");
  file << "H_J = [";
  for (j=0; j<total_steps; j++) file << H_j[j] << ", ";

  file << "]\nH_K = [";
  for (j=0; j<total_steps; j++) file << H_k[j] << ", ";
  file << "]";


}




















//
//
//
// #include <fstream>
// #include <iostream>
// #include <stdio.h>
// #include <stdlib.h>
// #include <string.h>
// #include <complex.h>
// #include <math.h>
// #include <gsl/gsl_rng.h>
// #include <string>
// #include <algorithm>
// #include <ctime>
//
// #include "check.h"
// #include "hamiltonian.h"
// #include "operations.h"
// #include "parameters.h"
// #include "print.h"
// #include "linear_algebra.h"
//
//
//
// main (int argc, char *argv[]){
//
//   Simulation_Parameters sim_params;
// 	sim_params.initialize_lattice(2,sim_params);
//
//
//   double tau = 0.351562, time_step = 0.0001, time = tau;
//   int i,j;
// 	double re =0, im = 0;
//   int total_steps = (int)floor(tau/time_step), N = 36;
//   double *H_j, *H_k, *pi, *psi, *psi_temp1, *psi_temp2, *ham_init, *ham_target, *target_state, *initial_state,*evolved_state;
//   double jkb_initial[3]     = {0.594,0.125, 0};
//   double jkb_target[3]     = {0.291,0.5,0};
//   initial_state           = new double[2*N]();
//   target_state           = new double[2*N]();
//   evolved_state           = new double[2*N]();
//   ham_init           = new double[2*N*N]();
//   ham_target           = new double[2*N*N]();
//
//   construct_device_hamiltonian_uniform(sim_params, ham_init, jkb_initial);
//   construct_device_hamiltonian_uniform(sim_params, ham_target, jkb_target);
// 	get_ground_state(N, ham_init, initial_state);
// 	get_ground_state(N, ham_target, target_state);
//   //
//   // double target_state[N*2] = {-0.0392708, 0, -0.0354592, 0, -0.0392708, 0, -0.136812, 0, -0.244756, 0, -0.0392708, 0, -0.244756, 0, -0.198235, 0, -0.244756, 0, -0.198235, 0, -0.0392708, 0, -0.244756, 0, -0.136812, 0, -0.308535, 0, -0.198235, 0, -0.0354592, 0, -0.136812, 0, -0.0551019, 0, -0.136812, 0, -0.244756, 0, -0.0392708, 0, -0.136812, 0, -0.308535, 0, -0.136812, 0, -0.244756, 0, -0.198235, 0, -0.244756, 0, -0.0392708, 0, -0.0551019, 0, -0.136812, 0, -0.0354592, 0, -0.0392708, 0, -0.244756, 0, -0.136812, 0, -0.0354592, 0, -0.0392708, 0};
//   // double evolved_state[N*2] = {0.0050459, 0.0358555, 0.0492864, 0.0336511, 0.0050459, 0.0358555, 0.0505893, 0.153499, 0.0647756, 0.227186, 0.0050459, 0.0358555, 0.0647756, 0.227186, 0.0234769, 0.189666, 0.0647756, 0.227186, 0.0234769, 0.189666, 0.0050459, 0.0358555, 0.0647756, 0.227186, 0.0505893, 0.153499, 0.0533961, 0.279351, 0.0234769, 0.189666, 0.0492864, 0.0336511, 0.0505893, 0.153499, 0.0449853, 0.062841, 0.0505893, 0.153499, 0.0647756, 0.227186, 0.0050459, 0.0358555, 0.0505893, 0.153499, 0.0533961, 0.279351, 0.0505893, 0.153499, 0.0647756, 0.227186, 0.0234769, 0.189666, 0.0647756, 0.227186, 0.0050459, 0.0358555, 0.0449853, 0.062841, 0.0505893, 0.153499, 0.0492864, 0.0336511, 0.0050459, 0.0358555, 0.0647756, 0.227186, 0.0505893, 0.153499, 0.0492864, 0.0336511, 0.0050459, 0.0358555};
//   H_j           = new double[total_steps]();
//   H_k           = new double[total_steps]();
//   pi           = new double[2*N]();
//   psi           = new double[2*N]();
//   psi_temp1      = new double[2*N]();
//   psi_temp2       = new double[2*N]();
//
//   print_state(initial_state, N);
//
//   for(i=0;i<N*2;i++){
//     psi[i] = evolved_state[i];
//   }
//
// 	for(i=0;i<N*2;i+=2){
// 		re += target_state[i]*evolved_state[i] + target_state[i+1]*evolved_state[i+1];
// 		im += target_state[i]*evolved_state[i+1] - target_state[i+1]*evolved_state[i];
// 	}
//   for(i=0;i<N*2;i+=2){
//     pi[i] = -2 * (re*target_state[i] - im*target_state[i+1]);
//     pi[i+1] = -2 * (re*target_state[i+1] + im*target_state[i]);
//
//   }
//
//
//
//   double *hamiltonian,*ham_t_i, *ham_real,*exp_matrix, *hamiltonian10, *hamiltonian01;
//
//   hamiltonian = new double[2*sim_params.N*sim_params.N]();
//   hamiltonian10 = new double[2*sim_params.N*sim_params.N]();
//   hamiltonian01 = new double[2*sim_params.N*sim_params.N]();
// 	exp_matrix  = new double[2*sim_params.N*sim_params.N]();
// 	ham_t_i     = new double[2*sim_params.N*sim_params.N]();
// 	ham_real    = new double[sim_params.N*sim_params.N]();
//
//
//
//
//
//









//
//
//
//   double jkb_vals[3] = {0,0,0};
//   double j_opt[1] = {0.00948312};
//   double k_opt[4] = {0, 0.0983751, 0.134202, 0.316601};
//   int j_index = 0;
//   int k_index = 0;
//   time = 0;
//   while(time <= 0.351562){
//     ///manage jkb_vals
//     if (j_index <= 1) {
//       if(time > j_opt[j_index]){
//         jkb_vals[0] = fmod(jkb_vals[0] +1, 2.0);
//         j_index += 1;
//         printf("time : %f, j: %f, k: %f\n", time, jkb_vals[0], jkb_vals[1]);
//       }
//     }
//     if (k_index <= 3) {
//       if(time > k_opt[k_index]){
//         jkb_vals[1] = fmod(jkb_vals[1] +1, 2.0);
//         k_index -= 1;
//         printf("time : %f, j: %f, k: %f\n", time, jkb_vals[0], jkb_vals[1]);
//       }
//     }
//     construct_device_hamiltonian_uniform(sim_params, hamiltonian, jkb_vals);
// 		for (j=0; j<sim_params.N*sim_params.N; j++) ham_t_i[2*j+1] = (hamiltonian[2*j]*time_step), ham_t_i[2*j] = -(hamiltonian[2*j+1]*time_step);  //multiplying by i*dt for the Pade approximation
// 		exp_complex_double_matrix_pade(sim_params.N, ham_t_i, exp_matrix);
// 		if(CHECK) check_unitary(exp_matrix, sim_params.N);
//
//     matrix_vector_mult(exp_matrix, psi, sim_params.N);
//
//
//     if(CHECK) check_norm(psi, N);
//
//
//     time  += time_step;
//   }
//
//   exit(0);
//
//   jkb_vals[3] = {1,0,0};
//   j_on     = {1,0,0};
//   k_on    = {0,1,0};
//   j_index = 0;
//   k_index = 3;
//   double j_on[3]     = {1,0,0};
//   double k_on[3]     = {0,1,0};
//   int h_index = total_steps-1;
//
//   construct_device_hamiltonian_uniform(sim_params, hamiltonian10, j_on);
//   construct_device_hamiltonian_uniform(sim_params, hamiltonian01, k_on);
//
//   while(time >= 0){
//     ///manage jkb_vals
//     if (j_index >= 0) {
//       if(time < j_opt[j_index]){
//         jkb_vals[0] = fmod(jkb_vals[0] +1, 2.0);
//         j_index -= 1;
//         //printf("time : %f, j: %f, k: %f\n", time, jkb_vals[0], jkb_vals[1]);
//       }
//     }
//     if (k_index >= 0) {
//       if(time < k_opt[k_index]){
//         jkb_vals[1] = fmod(jkb_vals[1] +1, 2.0);
//         k_index -= 1;
//         //printf("time : %f, j: %f, k: %f\n", time, jkb_vals[0], jkb_vals[1]);
//       }
//     }
//     construct_device_hamiltonian_uniform(sim_params, hamiltonian, jkb_vals);
// 		for (j=0; j<sim_params.N*sim_params.N; j++) ham_t_i[2*j+1] = (hamiltonian[2*j]*time_step), ham_t_i[2*j] = -(hamiltonian[2*j+1]*time_step);  //multiplying by i*dt for the Pade approximation
// 		exp_complex_double_matrix_pade(sim_params.N, ham_t_i, exp_matrix);
// 		if(CHECK) check_unitary(exp_matrix, sim_params.N);
//
//     matrix_vector_mult(exp_matrix, pi, sim_params.N);
//     matrix_vector_mult(exp_matrix, psi, sim_params.N);
//
//
//     memcpy(psi_temp1, psi, sizeof(double)*2*N);
//     memcpy(psi_temp2, psi, sizeof(double)*2*N);
//
//     if(CHECK) check_norm(psi, N);
//
//     construct_device_hamiltonian_uniform(sim_params, hamiltonian10, j_on);
//     construct_device_hamiltonian_uniform(sim_params, hamiltonian01, k_on);
//
//     matrix_vector_mult(hamiltonian10, psi_temp1, sim_params.N);
//     matrix_vector_mult(hamiltonian01, psi_temp2, sim_params.N);
//
//
//     for (i=0; i<2*N; i+=2){
//       H_j[h_index] += pi[i] * (psi_temp1[i+1])   - pi[i+1] * (psi_temp1[i]);
//       H_k[h_index] += pi[i] * (psi_temp2[i+1])   - pi[i+1] * (psi_temp2[i]);
//     }
//     h_index -= 1;
//     time  -= time_step;
//   }
//   print_state(psi, N);
//
//
//
//
//   std::ofstream file;
// 	file.open("../data_paper_98/h_data.txt");
//   file << "H_J = [";
//   for (j=0; j<total_steps; j++) file << H_j[j] << ", ";
//
//   file << "]\nH_K = [";
//   for (j=0; j<total_steps; j++) file << H_k[j] << ", ";
//   file << "]";
//
//
// }
