#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <complex>
#include <math.h>
#include <sys/time.h>
#include <gsl/gsl_rng.h>
#include <sys/stat.h>
#include <string>
#include <algorithm>

//g++ -o qubit_simulation qubit_simulation.cpp -llapack -lblas -lgsl
//./qubit_simulation

//Physical Parameters
const bool PERIODIC = true;
const bool UNIFORM_SITES = true;
const int DEVICE_DIMENSION = 2;
const double MAX_PARAM = 1;
const double MIN_PARAM = 0;
const double UPJ = MAX_PARAM;
const double UPK = MAX_PARAM;
const double UPB = MAX_PARAM;
const double LOWJ = MIN_PARAM;
const double LOWK = MIN_PARAM;
const double LOWB = MIN_PARAM;
const double T = 1;
const double V = 2;
const int NX = 3;
const int NY = 3;
const int NUM_SITES = NX*NY;


//Non-Physical Parameters
const bool MC = false;
const bool MC_DATA = false;
const bool MC_BANG = true;
const bool MCBB_DATA = true;
const bool ADIABATIC_DATA = false;
const bool ADIABATIC = false;
const bool EVOLVE_DATA = false;


//Parameters for both simulations
const bool DIAG = false;
const int SEED_TOTAL= 8;


//MC method parameters
const double DIFFERENCE_LIMIT_MC = 0.005;
const double TAU_INIT_MC = 0.03;
const double MAX_TAU_MC = 10;
const double TAU_SCALAR_MC = 1.4;
const double MAX_CHANGE_MC_INIT = 0.2*(MAX_PARAM-MIN_PARAM);
const double ACCEPTANCE_PROB_MC = 0.8;
const double TEMP_EXP_DECAY_MC = 0.85;
const double BINARY_SEARCH_TAU_LIMIT_MC = TAU_INIT_MC/30.0;
const int RANDOM_STATES_MC = 5;
const int SWEEPS_MC = 50;
const int TOTAL_STEPS_INIT_MC =  6;
const int TEMP_DECAY_ITERATIONS_MC =  30;
const int TEMP_DECAY_LIMIT_MC = 10;
const int MAX_EVOLVE_STEPS_MC = ceil((MAX_TAU_MC*TOTAL_STEPS_INIT_MC)/TAU_INIT_MC);
const int MAX_TAU_STEPS_MC = ceil(log(MAX_TAU_MC/TAU_INIT_MC)/log(TAU_SCALAR_MC));
const int ARRAY_SCALAR = 2;
double CONVERGE_LIMITS [6] = {.9,.92,.94,.96,.98, 1};


//MC_bang method parameters
const double DIFFERENCE_LIMIT_MCBB = DIFFERENCE_LIMIT_MC;
const double TAU_INIT_MCBB = TAU_INIT_MC;
const double MAX_TAU_MCBB = MAX_TAU_MC;
const double TAU_SCALAR_MCBB = TAU_SCALAR_MC;
const double CHANGE_FRACTION_MCBB = 0.1;
const double ACCEPTANCE_PROB_MCBB = ACCEPTANCE_PROB_MC;
const double TEMP_EXP_DECAY_MCBB = TEMP_EXP_DECAY_MC;
const double BINARY_SEARCH_TAU_LIMIT_MCBB = TAU_INIT_MCBB/20.0;
const int RANDOM_STATES_MCBB = RANDOM_STATES_MC;
const int NUMBER_OF_BANGS = 5;
const int SWEEPS_MCBB = 50;
const int TEMP_DECAY_ITERATIONS_MCBB = TEMP_DECAY_ITERATIONS_MC;;
const int TEMP_DECAY_LIMIT_MCBB = TEMP_DECAY_LIMIT_MC;
const int MAX_TAU_STEPS_MCBB = ceil(log(MAX_TAU_MCBB/TAU_INIT_MCBB)/log(TAU_SCALAR_MCBB));


//Adiabatic method parameters
const double DIFFERENCE_LIMIT_ADIA = 0.001;
const double TAU_INIT_ADIA = 0.15;
const double MAX_TAU_ADIA = 10;
const double TAU_SCALAR_ADIA = 1.3;
const double TIME_STEP_ADIA = 1/10000.0;





//Debugging Help
const bool CHECK = true;
const bool PRINT = true;
const bool PRINT_MC_ITERATIONS = true;
const bool PRINT_COMMUTATOR = false;
const bool PRINTBEST = false;
const bool PRINT_TIMES = false;

//The change_array variables. If 0 -> This method will not be used, if n -> use this method on every nth iteration of change_array
const int ROW = 0;
const int COL = 0;
const int ALTROW = 0;
const int ALTCOL = 0;
const int ALTROW2 = 0;
const int ALTCOL2 = 0;
const int ROWK =1;
const int ROWJ =2;
const int ROWB =3;
const int COLK =0;
const int COLJ =0;
const int COLB =0;
const int SINGLE =0;
const int VARIATIONS = (ROW && 1)+(COL && 1)+(ALTROW && 1)+(ALTCOL && 1)+(ALTROW2 && 1)+(ALTCOL2 &&  1)+(SINGLE &&  1)+(ROWK && 1)+(ROWJ && 1)+(ROWB && 1)+(COLK && 1)+(COLJ && 1)+(COLB && 1);




using namespace std;


//Homemade functions
void monte_carlo_method(gsl_rng * r,int *table, unsigned long long int *b, int num_electrons, int N, int lattice[NX][NY], double*ham_target,double* psi_start,double* jkb_initial,double* jkb_target,double ground_E, double *gf);

int ground_state_binary_search(int *table, unsigned long long int *b, int num_electrons, int N, int lattice[NX][NY], double*ham_target,int total_steps, double time_step,double* psi_start,gsl_rng * r,int seed,double* jkb_initial, double* jkb_target,double ground_E, double tau, double tau_min_bs,double* k_best,double* j_best,double* b_best,double* best_E_array,double* tau_array,int index);

int ground_state_binary_search_bang_bang(int *table, unsigned long long int *b, int num_electrons, int N, int lattice[NX][NY], double*ham_target, double* psi_start,gsl_rng * r,int seed,double* jkb_initial, double* jkb_target,double ground_E, double tau, double tau_min_bs,double* j_best_times,double* k_best_times,double* b_best_times,double* best_E_array,double* tau_array,int index, double *pre_converge_Es, double *pre_converge_taus, double *pre_converge_j_times,double * pre_converge_k_times,double *pre_converge_b_times);

double monte_carlo_simulation_fixed_tau(int *table, unsigned long long int *b, int num_electrons, int N, int lattice[NX][NY], double*ham_target, double* psi_start, double temp, double tau, int total_steps, double time_step, gsl_rng * r, double *k_best,double *j_best, double *b_best, double ground_E, double *best_E_array, int seed);

void evolve_mc(int *table, unsigned long long int *b,int num_electrons,int N,  int lattice[NX][NY], double *psi, double *k_array, double *j_array, double *b_array, int *bonds, int total_steps, double time_step);

void monte_carlo_method_bang_bang(gsl_rng * r,int *table, unsigned long long int *b, int num_electrons, int N, int lattice[NX][NY], double*ham_target,double* psi_start,double* jkb_initial,double* jkb_target,double ground_E, double *gf);

void check_pre_convergence(double *pre_converge_Es,double *pre_converge_taus,double *pre_converge_j_times, double *pre_converge_k_times, double *pre_converge_b_times, double *j_best_times,double *k_best_times,double *b_best_times,double best_E, double ground_E, double tau);

double monte_carlo_simulation_bang_bang_fixed_tau(int *table, unsigned long long int *b, int num_electrons, int N, int lattice[NX][NY], double*ham_target, double* psi_start, double temp, double tau, gsl_rng * r, double *j_best_times,double *k_best_times, double *b_best_times, double ground_E, double *best_E_array, int seed);

void evolve_mc_bang_bang(int *table, unsigned long long int *b,int num_electrons,int N,  int lattice[NX][NY], double *psi, double *j_times, double *k_times, double *b_times, double tau);

void print_times(double* j_times, double* k_times, double* b_times);

void adiabatic_method(gsl_rng * r,int *table, unsigned long long int *b, int num_electrons, int N, int lattice[NX][NY], double*ham_target,double* psi_start,double* jkb_initial,double* jkb_target,double g_initial, double f_initial, double g_target, double f_target,double ground_E);

void evolve_adiabatic(int *table, unsigned long long int *b,int num_electrons,int N, int lattice[NX][NY], double *psi, int total_steps, double time_step,double tau, double g_i, double f_i, double g_t, double f_t);

double cost(double *psi, double *ham_target, int N);

double calc_initial_temp_bang_bang(int *table, unsigned long long int *b, int num_electrons, int N, int lattice[NX][NY], double*ham_target, double tau, double *psi_start, gsl_rng *r, int seed, double *jkb_initial, double *jkb_target);

double calc_initial_temp(int *table, unsigned long long int *b, int num_electrons, int N, int lattice[NX][NY], double*ham_target, int total_steps, double time_step, double tau, double *psi_start, gsl_rng *r, int seed, double *jkb_initial, double *jkb_target);

double get_change_mc(gsl_rng * r, double tau);

double get_change_mcbb(gsl_rng * r, double tau);

void construct_device_hamiltonian(int *table, unsigned long long int *b,int num_electrons,int N, double *ham_dev, int lattice[NX][NY],double *k_array, double *j_array, double *b_array, int *bonds, int index, int D);

void construct_device_hamiltonian_uniform(int *table, unsigned long long int *b,int num_electrons,int N, double *ham_dev, int lattice[NX][NY],double *jkb, int D);

void construct_model_hamiltonian(int *table, unsigned long long int *b,int num_electrons,int N, double *ham_target,  int lattice[NX][NY], int D);

double find_next_time(double time, double tau,double* j_times, double* k_times,double*  b_times,int*  jkb_index, double* jkb);

void diag_hermitian_real_double(int N,  double *A, double *Vdag,double *D);

void exp_diaganolized_mat(double *ham_real, double *Vdag, double* D, int N, double time_step);

void exp_general_complex_double(int N, double *A, double *B);

void matrix_vector_mult(double *exp_matrix, double *psi, int N);

int hop(unsigned long long int b, unsigned long long int *v,int n, int j);

int find(int N,unsigned long long int *v, unsigned long long int *b);

void get_ground_state(int N, double *ham, double* ground_state);

double get_ground_E(int N, double *ham);

unsigned long long int choose(int num_electrons);

int combinations ( int num_electrons, unsigned long long int *b,int *tab, int N);

void assign_bonds(int *bonds, int lattice[NX][NY]);

void copy_arrays(int N, double *k_array, double *j_array, double* b_array,  double* k_best,  double* j_best, double* b_best, int total_steps);

void copy_arrays_bang_bang(double *j_to, double *k_to, double* b_to,  double* j_from,  double* k_from, double* b_from);

void init_arrays(double *k_array, double *j_array,double *b_array, gsl_rng *r, int total_steps);

void init_arrays_bang(double *j_best_times, double *k_best_times,double *b_best_times, gsl_rng *r);

void change_array(double *k_array, double *j_array, double *b_array, int random_row, int random_col, double change, int i, int total_steps);

void change_array_bang_bang(double *j_times, double *k_times, double *b_times, double change, int random_time_index, int i, double tau);

void change_row(double *k_array, double *j_array,double *b_array, int row, double change, bool k, bool j, bool b, int jump, int offset);

void change_col(int total_steps,double *k_array,double *j_array,double *b_array, int col, double change,bool k, bool j, bool b, int jump, int offset);

void change_single(double *k_array,double *j_array,double *b_array, int row,int col, double change);

void scale_arrays_bang_bang(double *j, double *k, double *b, double scalar);

void scale_arrays(int total_steps, double* j_best, double* k_best,double* b_best);

void construct_lattice(int lattice[NX][NY]);

int get_neighbors(int site, int *neighbors, int lattice[NX][NY], int D);

double get_random(double lower, double upper, gsl_rng *r);

void check_commutator(int N, double* A, double* B);

void check_norm(double* psi, int N);

void check_unitary(double* ham, int N);

void check_hermicity(double* ham_target, int N);

void check_weights(double* state, double* ham, int N);

void export_mc_data(double pre_converge_tau, int pre_converge_steps, double pre_converge_E, double *j_best, double *k_best, double *b_best,double *pre_converge_j_best, double *pre_converge_k_best, double *pre_converge_b_best, double *tau_array, double *best_E_array, int total_steps, double tau, double ground_E, double *jkb_initial, double *jkb_target, int index, int seed, int num_electrons, double *gf);

void export_mcbb_data(double* pre_converge_taus,double* pre_converge_Es,double* pre_converge_j_best,double*  pre_converge_k_best,double* pre_converge_b_best, double *j_best_times, double *k_best_times, double *b_best_times, double *tau_array, double *best_E_array, double tau, double ground_E, double *jkb_initial, double *jkb_target, int index, int seed, int num_electrons, double *gf);

void export_adiabatic_data(double *tau_array, double *E_array, int total_steps, double tau, double ground_E, double g_initial, double f_initial, double g_target, double f_target,double *jkb_initial, double *jkb_taget, int index, int num_electrons);

void export_evolve_data(int *table, unsigned long long int *b, int num_electrons, int N,double* j_best, double* k_best, double* b_best, double tau, double time_step, int total_steps, int lattice[NX][NY],double* psi_start, double ground_E, double* ham_target, int seed);

void export_evolve_calc(int *table, unsigned long long int *b,int num_electrons,int N,  int lattice[NX][NY], double *psi, double *k_array, double *j_array, double *b_array, int *bonds, int total_steps, double time_step, double* E_list, double* T_list, double* ham_target);

void print_best_arrays(double *k_array, double *j_array, double *b_array, int totals_step);

void print_vector(double* psi,int N);

void print_hamiltonian(double* hamiltonian, int N);

void print_hamiltonianR(double *hamiltonian, int N);

void print_E(double ground_E, double best_E, double initial_E);

void test_function();


//All of the linear algebra functions, courtesy of llapack/llblas
extern "C" int zgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K, double *ALPHA, double *Z, int *LDA, double *X, int *LDB, double *BETA, double *Y, int *LDC); //complex matrix*matrix mult, odd indices hold imaginary values.
extern "C" int zgemv_(char *TRANS, int *M, int *N,double *ALPHA,double *A, int *LDA, double *X, int *INCX, double *BETA, double *Y, int *INCY); //complex matrix-vector mult, odd indices hold imaginary values.
extern "C" int dsyev_(char *JOBZ, char *UPLO, int *N, double *Vdag, int *LDA, double *D, double *WORK, int *LWORK, int *INFO);//diagonalization, returns the eigenvectors in Vdag and eigenvalues in D.
extern "C" int zgetrf_ (int *M, int *N, double *D, int *LDA, int *IPIV, int *INFO);//A matrix factorization?
extern "C" int zgetrs_( char *TRANS, int *N, int *NRHS,double *D, int *LDA, int *IPIV,double *E, int *LDB, int *INFO);//solves a system of equations
extern "C" void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
extern "C" void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
extern "C" int dgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K, double *ALPHA, double *Z, int *LDA, double *X, int *LDB, double *BETA, double *Y, int *LDC); //real values matrix*matrix mult
extern "C" double zdotc_(int *N, double*ZX,int *INCX, double *ZY, int *INCY);//dots the complex conjugate of ZX with ZY









main (int argc, char *argv[])
{
	int *table,lattice[NX][NY],num_electrons,N, gt_i, ft_i, gi_i, fi_i, electron_index;
	unsigned long long int *b;
	double *ham_target, *ham_initial, *psi_start, ground_E, *jkb_initial, *jkb_target, g_initial, f_initial, g_target, f_target;

	construct_lattice(lattice);

	gsl_rng_env_setup();
	const gsl_rng_type * TT = gsl_rng_default;
	gsl_rng * r  = gsl_rng_alloc (TT);
	jkb_initial = (double *) malloc(3*sizeof(double));
	jkb_target = (double *) malloc(3*sizeof(double));


	double f_targ[3] = {0.5,0.1,1};
	double g_targ[3] = {0.5,0.1,1};
	double f_init[2] = {0.1,1};
	double g_init[2] = {0.1,1};

	int num_electrons_array[4]= {2,3,4,7};
	for(electron_index=0;electron_index<sizeof(num_electrons_array)/sizeof(int);electron_index++)
	{
		num_electrons = num_electrons_array[electron_index];
		N=choose(num_electrons);
		b = (unsigned long long int*) malloc(N*N*N*sizeof(unsigned long long int));
		table=(int*) malloc(num_electrons*N*sizeof(int));
		ham_target = (double *) malloc (2*N*N*sizeof (double));
		ham_initial = (double *) malloc(2*N*N*sizeof(double));
		psi_start = (double *) malloc(2*N*sizeof(double));
		combinations (num_electrons,b,table, N);

		for(ft_i = 0; ft_i < (sizeof(f_targ) / sizeof(double)); ft_i++)
		{
			f_target = f_targ[ft_i];
			for(gt_i=0; gt_i < (sizeof(g_targ) / sizeof(double)) ; gt_i++)
			{
				g_target = g_targ[gt_i];
				for(fi_i = 0; fi_i < (sizeof(f_init) / sizeof(double)) ; fi_i++)
				{
					f_initial = f_init[fi_i];
					for(gi_i = 0; gi_i < (sizeof(g_init) / sizeof(double)) ; gi_i++)
					{
						g_initial = g_init[gi_i];

						jkb_initial[0] = g_initial;
						jkb_initial[1] = (1-g_initial)*f_initial;
						jkb_initial[2] = (1-f_initial)*(1-g_initial);
						construct_device_hamiltonian_uniform(table, b, num_electrons, N, ham_initial, lattice, jkb_initial, DEVICE_DIMENSION);
						printf("\n\nHAM_INITIAL:\n");
						print_hamiltonian(ham_initial, N);

						jkb_target[0] = g_target;
						jkb_target[1] = (1-g_target)*f_target;
						jkb_target[2] = (1-f_target)*(1-g_target);
						construct_device_hamiltonian_uniform(table, b, num_electrons, N, ham_target, lattice, jkb_target, DEVICE_DIMENSION);
						printf("\n\nHAM_TARGET:\n");
						print_hamiltonian(ham_target, N);

						if(CHECK) check_commutator(N, ham_initial, ham_target);

						get_ground_state(N, ham_initial,psi_start);
						ground_E = get_ground_E(N, ham_target);
						if(CHECK) check_norm(psi_start, N);




						double gf[4] = {g_initial, f_initial, g_target, f_target};


						if(MC) monte_carlo_method(r,table, b, num_electrons, N,lattice,ham_target, psi_start, jkb_initial,jkb_target, ground_E, gf);
						if(MC_BANG) monte_carlo_method_bang_bang(r,table, b, num_electrons, N,lattice,ham_target, psi_start, jkb_initial,jkb_target, ground_E, gf);
						if(ADIABATIC)  adiabatic_method(r,table, b, num_electrons, N,lattice,ham_target, psi_start, jkb_initial,jkb_target,g_initial, f_initial, g_target, f_target, ground_E);
					}//g_init closed bracket
				}//f_init closed bracket
			}//g_target close bracket
		}//f_target close bracket
		free(ham_initial), free(ham_target), free(b), free(table), free(psi_start);
	}
}








void monte_carlo_method(gsl_rng * r,int *table, unsigned long long int *b, int num_electrons, int N, int lattice[NX][NY], double*ham_target,double* psi_start,double* jkb_initial,double* jkb_target,double ground_E, double *gf)
{
	int seed, i, index, total_steps, pre_converge_steps;
	double initial_temp, *k_best, *j_best, *b_best, *pre_converge_k_best, *pre_converge_j_best, *pre_converge_b_best,*best_E_array, *tau_array, difference, tau, tau_min_bs, time_step, pre_converge_tau, pre_converge_E, best_E;

	k_best = (double *) malloc(2*NUM_SITES*MAX_EVOLVE_STEPS_MC*sizeof(double));
	j_best = (double *) malloc(2*NUM_SITES*MAX_EVOLVE_STEPS_MC*sizeof(double));
	b_best = (double *) malloc(NUM_SITES*MAX_EVOLVE_STEPS_MC*sizeof(double));
	pre_converge_k_best = (double *) malloc(2*NUM_SITES*MAX_EVOLVE_STEPS_MC*sizeof(double));
	pre_converge_j_best = (double *) malloc(2*NUM_SITES*MAX_EVOLVE_STEPS_MC*sizeof(double));
	pre_converge_b_best = (double *) malloc(NUM_SITES*MAX_EVOLVE_STEPS_MC*sizeof(double));
	tau_array = (double *) malloc(MAX_TAU_STEPS_MC*sizeof(double));
	best_E_array = (double *) malloc(MAX_TAU_STEPS_MC*sizeof(double));

	for(seed=1; seed<SEED_TOTAL+1; seed++)
	{
		gsl_rng_set(r, seed);
		for (i=0;i<2*NUM_SITES*MAX_EVOLVE_STEPS_MC;i++) k_best[i] =0, j_best[i] = 0;
		for (i=0;i<NUM_SITES*MAX_EVOLVE_STEPS_MC;i++) b_best[i] =0;

		index = 0;
		init_arrays(k_best, j_best, b_best, r, MAX_EVOLVE_STEPS_MC);

		tau = TAU_INIT_MC, tau_min_bs = TAU_INIT_MC;
		total_steps = TOTAL_STEPS_INIT_MC;
		time_step = tau/((double) total_steps);
		pre_converge_tau = 0,  pre_converge_E = 0, pre_converge_steps = 0;

		while(tau<MAX_TAU_MC)
		{
			initial_temp = calc_initial_temp(table, b, num_electrons, N, lattice, ham_target,total_steps, time_step,tau, psi_start, r, seed, jkb_initial, jkb_target);
			best_E = monte_carlo_simulation_fixed_tau(table, b, num_electrons, N, lattice, ham_target, psi_start, initial_temp,tau, total_steps, time_step,r, k_best, j_best, b_best, ground_E, best_E_array, seed);

			difference = best_E-ground_E;

			if((difference < (best_E_array[0]-ground_E)/10.0) && pre_converge_tau==0) pre_converge_tau=tau, pre_converge_E = best_E, pre_converge_steps = total_steps, copy_arrays(N, k_best, j_best,b_best, pre_converge_k_best, pre_converge_j_best, pre_converge_b_best, total_steps);


			if(difference < DIFFERENCE_LIMIT_MC)
			{
				index = ground_state_binary_search(table, b, num_electrons, N, lattice, ham_target,total_steps, time_step, psi_start, r, seed, jkb_initial, jkb_target, ground_E, tau, tau_min_bs, k_best, j_best, b_best, best_E_array, tau_array, index);
				break;
			}
			else
			{
				tau_array[index] = tau;
				best_E_array[index] = best_E;
				index ++;
				tau_min_bs = tau;
				tau = tau*TAU_SCALAR_MC;
				total_steps = floor(TOTAL_STEPS_INIT_MC*(tau/TAU_INIT_MC));
				time_step = tau/((double) total_steps);

			}
		}
		if(MC_DATA) export_mc_data(pre_converge_tau, pre_converge_steps, pre_converge_E, j_best, k_best,b_best,pre_converge_k_best, pre_converge_j_best, pre_converge_b_best,tau_array, best_E_array, total_steps,tau, ground_E, jkb_initial, jkb_target, index, seed, num_electrons, gf);
	}
	free(k_best),free(j_best), free(b_best), free(pre_converge_k_best), free(pre_converge_j_best), free(pre_converge_b_best), free(tau_array), free(best_E_array);
}







int ground_state_binary_search(int *table, unsigned long long int *b, int num_electrons, int N, int lattice[NX][NY], double*ham_target,int total_steps, double time_step,double* psi_start,gsl_rng * r,int seed,double* jkb_initial, double* jkb_target,double ground_E, double tau, double tau_min_bs,double* k_best,double* j_best,double* b_best,double* best_E_array,double* tau_array,int index)
{
	double tau_max_bs, initial_temp, best_E, difference;
	tau_max_bs = tau;
	tau = (tau_max_bs + tau_min_bs) / 2.0;

	printf("\nUsing binary search method to look for optimal ground state...");

	while((tau_max_bs - tau_min_bs) >BINARY_SEARCH_TAU_LIMIT_MC)
	{
		initial_temp = calc_initial_temp(table, b, num_electrons, N, lattice, ham_target,total_steps, time_step, tau, psi_start, r, seed, jkb_initial, jkb_target);
		best_E = monte_carlo_simulation_fixed_tau(table, b, num_electrons, N, lattice, ham_target, psi_start, initial_temp,tau, total_steps, time_step,r, k_best, j_best, b_best, ground_E, best_E_array, seed);
		difference = best_E-ground_E;

		if(difference < DIFFERENCE_LIMIT_MC)//still in the ground state, backtrack tau
		{
			printf("\nIN THE BINARY SEARCH PROCESS\nSTEPPING BACKWARD....\n");
			tau_max_bs = tau;
			tau = (tau_max_bs+tau_min_bs)/2.0;
			time_step = tau/((double) total_steps);
		}
		else
		{
			printf("\nIN THE BINARY SEARCH PROCESS\nSTEPPING FORWARD...\n");
			tau_array[index] = tau;
			best_E_array[index] = best_E;
			index++;
			tau_min_bs = tau;
			tau = (tau_min_bs+tau_max_bs)/2;
			time_step = tau/((double) total_steps);
		}
	}
	tau_array[index] = tau;
	best_E_array[index] = best_E;
	printf("\nGROUND STATE REACHED, EXITING BINARY SEARCH PROCESS\n");
	return index;
}


int ground_state_binary_search_bang_bang(int *table, unsigned long long int *b, int num_electrons, int N, int lattice[NX][NY], double*ham_target, double* psi_start,gsl_rng * r,int seed,double* jkb_initial, double* jkb_target,double ground_E, double tau, double tau_min_bs,double* j_best_times,double* k_best_times,double* b_best_times,double* best_E_array,double* tau_array,int index, double *pre_converge_Es, double *pre_converge_taus, double *pre_converge_j_times,double * pre_converge_k_times,double *pre_converge_b_times)
{
	double tau_max_bs, initial_temp, best_E, difference;
	bool step_backward = true;

	tau_max_bs = tau;
	tau = (tau_max_bs + tau_min_bs) / 2.0;


	printf("\nUsing binary search method to look for optimal ground state...");

	while((tau_max_bs - tau_min_bs) >BINARY_SEARCH_TAU_LIMIT_MCBB)
	{
		if(step_backward) scale_arrays_bang_bang(j_best_times,k_best_times,b_best_times,tau/tau_max_bs);
		else scale_arrays_bang_bang(j_best_times,k_best_times,b_best_times,tau/tau_min_bs);

		initial_temp = calc_initial_temp_bang_bang(table, b, num_electrons, N, lattice, ham_target,tau, psi_start, r, seed, jkb_initial, jkb_target);
		best_E = monte_carlo_simulation_bang_bang_fixed_tau(table, b, num_electrons, N, lattice, ham_target, psi_start, initial_temp,tau, r, j_best_times, k_best_times, b_best_times, ground_E, best_E_array, seed);

		check_pre_convergence(pre_converge_Es, pre_converge_taus, pre_converge_j_times, pre_converge_k_times, pre_converge_b_times, j_best_times,k_best_times,b_best_times,best_E, ground_E,tau);

		difference = best_E-ground_E;

		if(difference < DIFFERENCE_LIMIT_MC)//still in the ground state, backtrack tau
		{
			printf("\nIN THE BINARY SEARCH PROCESS\nSTEPPING BACKWARD....\n");
			tau_max_bs = tau;
			tau = (tau_max_bs+tau_min_bs)/2.0;
			step_backward = true;
		}
		else
		{
			printf("\nIN THE BINARY SEARCH PROCESS\nSTEPPING FORWARD...\n");
			tau_array[index] = tau;
			best_E_array[index] = best_E;
			index++;
			tau_min_bs = tau;
			tau = (tau_min_bs+tau_max_bs)/2;
			step_backward = false;
		}
	}

	if(step_backward)tau_array[index] = tau_max_bs, best_E_array[index] = best_E;
	else index--;
	printf("\nGROUND STATE REACHED, EXITING BINARY SEARCH PROCESS\n");
	return index;
}




void scale_arrays_bang_bang(double *j, double *k, double *b, double scalar)
{
	int i;
	for(i=0;i<NUMBER_OF_BANGS*2;i++)
	{
		j[i] = j[i]*scalar;
		k[i] = k[i]*scalar;
		b[i] = b[i]*scalar;
	}
}

double monte_carlo_simulation_fixed_tau(int *table, unsigned long long int *b, int num_electrons, int N, int lattice[NX][NY], double*ham_target, double* psi_start, double temperature, double tau, int total_steps, double time_step, gsl_rng * r, double *k_best,double *j_best, double *b_best, double ground_E, double *best_E_array, int seed)
/*monte_carlo_simulation_fixed_tau the values in the j, b, and k list in order to produce the lowest energy (expectation value) between the final state (psi), produced by evolve, and the model hamiltonian. This is done by randomly selecting one row of each list, making a slight change, then determining if the new energy is lower than the old. If the new energy is greater than the old, keep with probability exp(delta_E/Temp)*/
{
	int i=0,j=0, random_row=0, random_col, proposal_accepted=0,proposal_count=0, poor_acceptance_count=0,*bonds;
	double *psi, *k_array, *j_array,*b_array,*k_temp, *j_temp, *b_temp, acceptance_rate=0, E_old=0, E_new=0,best_E=0, change=0;
	bonds = (int*) malloc(NUM_SITES*NUM_SITES*sizeof(int));
	k_array = (double *) malloc(2*NUM_SITES*total_steps*sizeof(double));
	j_array = (double *) malloc(2*NUM_SITES*total_steps*sizeof(double));
	b_array = (double *) malloc(NUM_SITES*total_steps*sizeof(double));
	k_temp = (double *) malloc(2*NUM_SITES*total_steps*sizeof(double));
	j_temp = (double *) malloc(2*NUM_SITES*total_steps*sizeof(double));
	b_temp = (double *) malloc(NUM_SITES*total_steps*sizeof(double));
	psi = (double *) malloc (2*N*sizeof(double));

	assign_bonds(bonds, lattice);
	memcpy(psi,psi_start, 2*N*sizeof(double));

	if(CHECK) check_norm(psi, N);
	if(PRINTBEST) printf("\nPrinting the best K-J-B arrays from Monte_Carlo"), print_best_arrays(k_best, j_best, b_best, total_steps);

	copy_arrays(N, k_best, j_best,b_best, k_array, j_array, b_array, total_steps);
	evolve_mc(table,b,num_electrons,N,lattice, psi, k_array,j_array,b_array,bonds, total_steps, time_step);
	best_E = cost(psi, ham_target, N);

	if(total_steps == TOTAL_STEPS_INIT_MC && tau==TAU_INIT_MC) best_E_array[0] = best_E;
	E_old = best_E;
	if(PRINT) printf("Pre-Monte_Carlo Expectation:   %f\n", best_E);

	for (i=0;i<TEMP_DECAY_ITERATIONS_MC;i++)
	{
		proposal_accepted = 0;
		proposal_count = 0;

		for (j=0; j<SWEEPS_MC*total_steps;j++)
		{
			copy_arrays(N, k_array, j_array, b_array, k_temp, j_temp,b_temp,total_steps);//a temporary array, used in the undoing of the changes

			change = get_change_mc(r, tau);
			random_row = floor(get_random(0,total_steps,r));
			random_col = floor(get_random(0,NUM_SITES*2,r));

			change_array(k_array,j_array,b_array,random_row,random_col,change,j, total_steps);
			memcpy(psi,psi_start, 2*N*sizeof(double));//resetting psi

			evolve_mc(table,b,num_electrons,N,lattice, psi, k_array,j_array,b_array, bonds, total_steps, time_step);
			E_new = cost(psi, ham_target, N);
			if (E_new<best_E) best_E=E_new, E_old=E_new,  copy_arrays(N, k_array, j_array, b_array, k_best, j_best,b_best, total_steps), poor_acceptance_count = 0;
			else if (E_new<=E_old) E_old=E_new;
			else if (get_random(0,1,r)<exp(-(E_new-E_old)/(temperature))) E_old=E_new, proposal_accepted++, proposal_count++;
			else copy_arrays(N, k_temp, j_temp, b_temp, k_array, j_array, b_array,total_steps),proposal_count++;//undoing the change

		}

		acceptance_rate = (double)proposal_accepted/proposal_count;

		if(PRINT_MC_ITERATIONS) printf("           Best Expectation:   %3.6f  ||  Acceptance Rate: %3.4f (%i/%i)\n",best_E,acceptance_rate,proposal_accepted, proposal_count);

		if(acceptance_rate<0.1) poor_acceptance_count++;
		else poor_acceptance_count = 0;

		if(poor_acceptance_count>TEMP_DECAY_LIMIT_MC)
		{
			printf("NO MC PROGRESS FOR %i TEMP MAX_CHANGE_MC_INIT ITERATIONS, TERMINATING\n", TEMP_DECAY_LIMIT_MC);
			break;
		}
		temperature=temperature*TEMP_EXP_DECAY_MC;
	}

	if (PRINTBEST) printf("\nPrinting the best K-J-B arrays from main"), print_best_arrays(k_best, j_best, b_best, total_steps);
	if(EVOLVE_DATA) export_evolve_data(table, b, num_electrons, N, j_best, k_best, b_best, tau, time_step, total_steps, lattice, psi_start, ground_E, ham_target, seed);

	free(k_array), free(j_array), free(b_array),free(k_temp), free(j_temp), free(b_temp),free(psi), free(bonds);
	return best_E;
}







void monte_carlo_method_bang_bang(gsl_rng * r,int *table, unsigned long long int *b, int num_electrons, int N, int lattice[NX][NY], double*ham_target,double* psi_start,double* jkb_initial,double* jkb_target,double ground_E, double *gf)
{
	int seed, i, index, total_steps;
	double initial_temp, *k_best_times, *j_best_times, *b_best_times, *pre_converge_k_times, *pre_converge_j_times, *pre_converge_b_times,*best_E_array, *tau_array, difference, tau, tau_min_bs, time_step, *pre_converge_taus, *pre_converge_Es, best_E;

	j_best_times = (double *) malloc(2*NUMBER_OF_BANGS*sizeof(double));
	k_best_times = (double *) malloc(2*NUMBER_OF_BANGS*sizeof(double));
	b_best_times = (double *) malloc(2*NUMBER_OF_BANGS*sizeof(double));
	tau_array = (double *) malloc(MAX_TAU_STEPS_MCBB*sizeof(double));
	best_E_array = (double *) malloc(MAX_TAU_STEPS_MCBB*sizeof(double));


	pre_converge_k_times = (double *) malloc((sizeof(CONVERGE_LIMITS))*2*NUMBER_OF_BANGS);
	pre_converge_j_times = (double *) malloc((sizeof(CONVERGE_LIMITS))*2*NUMBER_OF_BANGS);
	pre_converge_b_times = (double *) malloc((sizeof(CONVERGE_LIMITS))*2*NUMBER_OF_BANGS);
	pre_converge_Es = (double *) malloc(sizeof(CONVERGE_LIMITS));
	pre_converge_taus = (double *) malloc(sizeof(CONVERGE_LIMITS));


	pre_converge_taus[0] = TAU_INIT_MCBB;

	for(seed=1; seed<SEED_TOTAL+1 ;seed++)
	{
		gsl_rng_set(r, seed);
		for (i=0;i<2*NUMBER_OF_BANGS;i++) j_best_times[i] = 0,k_best_times[i] = 0,b_best_times[i] = 0;
		for (i=0;i<2*NUMBER_OF_BANGS;i++) best_E_array[i] = 0, tau_array[i] = 0;
		for (i=0;i<sizeof(CONVERGE_LIMITS)/sizeof(double);i++) pre_converge_taus[i] = 0, pre_converge_Es[i] = 0;
		for (i=0;i<((sizeof(CONVERGE_LIMITS))*2*NUMBER_OF_BANGS)/sizeof(double);i++) pre_converge_k_times[i] = 0, pre_converge_j_times[i] = 0, pre_converge_b_times[i] = 0;



		index = 0;
		init_arrays_bang(j_best_times, k_best_times, b_best_times, r);
		copy_arrays_bang_bang(pre_converge_j_times, pre_converge_k_times, pre_converge_b_times,j_best_times, k_best_times, b_best_times);

		tau = TAU_INIT_MCBB;
		tau_min_bs = TAU_INIT_MCBB;




		while(tau<MAX_TAU_MCBB)
		{
			initial_temp = calc_initial_temp_bang_bang(table, b, num_electrons, N, lattice, ham_target,tau, psi_start, r, seed, jkb_initial, jkb_target);

			best_E = monte_carlo_simulation_bang_bang_fixed_tau(table, b, num_electrons, N, lattice, ham_target, psi_start, initial_temp,tau, r, j_best_times, k_best_times, b_best_times, ground_E, best_E_array, seed);

			if(tau==TAU_INIT_MCBB) pre_converge_Es[0] = best_E_array[0];
			if(PRINT) print_E(ground_E, best_E, pre_converge_Es[0]);


			difference = best_E-ground_E;

			check_pre_convergence(pre_converge_Es, pre_converge_taus, pre_converge_j_times, pre_converge_k_times, pre_converge_b_times, j_best_times,k_best_times,b_best_times,best_E, ground_E,tau);



			if(difference < DIFFERENCE_LIMIT_MCBB)
			{
				index = ground_state_binary_search_bang_bang(table, b, num_electrons, N, lattice, ham_target, psi_start, r, seed, jkb_initial, jkb_target, ground_E, tau, tau_min_bs, j_best_times, k_best_times, b_best_times, best_E_array, tau_array, index, pre_converge_Es, pre_converge_taus, pre_converge_j_times, pre_converge_k_times,pre_converge_b_times);
				if(index == 0)	best_E_array[index] = best_E;
				tau = tau_array[index];
				break;
			}
			else
			{
				tau_array[index] = tau;
				best_E_array[index] = best_E;
				index ++;
				tau_min_bs = tau;
				tau = tau*TAU_SCALAR_MCBB;

				scale_arrays_bang_bang(j_best_times,k_best_times,b_best_times,TAU_SCALAR_MCBB);
			}
		}
		if(MCBB_DATA) export_mcbb_data(pre_converge_taus, pre_converge_Es, pre_converge_j_times, pre_converge_k_times, pre_converge_b_times, j_best_times, k_best_times,b_best_times,tau_array, best_E_array, tau, ground_E, jkb_initial, jkb_target, index, seed, num_electrons, gf);

	}
	free(k_best_times),free(j_best_times), free(b_best_times),free(pre_converge_k_times), free(pre_converge_j_times), free(pre_converge_b_times),free(tau_array), free(best_E_array), free(pre_converge_Es), free(pre_converge_taus);
}



void check_pre_convergence(double *pre_converge_Es,double *pre_converge_taus,double *pre_converge_j_times, double *pre_converge_k_times, double *pre_converge_b_times, double *j_best_times,double *k_best_times,double *b_best_times,double best_E, double ground_E, double tau)
{
	int i;
	double initial_E = pre_converge_Es[0];
	double ratio = 1 - (best_E-ground_E)/(initial_E-ground_E);
	for(i=0; i<5;i++)
	{
		if(ratio >= CONVERGE_LIMITS[i] && ratio <= CONVERGE_LIMITS[i+1])
		{
			if(pre_converge_taus[i+1] == 0 || (tau < pre_converge_taus[i+1]))
			{
				pre_converge_taus[i+1] = tau;
				pre_converge_Es[i+1] = best_E;

				memcpy(pre_converge_j_times + (2*NUMBER_OF_BANGS*(i+1)), j_best_times, 2*NUMBER_OF_BANGS*sizeof(double));
				memcpy(pre_converge_k_times + (2*NUMBER_OF_BANGS*(i+1)), k_best_times, 2*NUMBER_OF_BANGS*sizeof(double));
				memcpy(pre_converge_b_times + (2*NUMBER_OF_BANGS*(i+1)), b_best_times, 2*NUMBER_OF_BANGS*sizeof(double));
			}
			break;
		}
	}

}



double monte_carlo_simulation_bang_bang_fixed_tau(int *table, unsigned long long int *b, int num_electrons, int N, int lattice[NX][NY], double*ham_target, double* psi_start, double temp, double tau, gsl_rng * r, double *j_best_times,double *k_best_times, double *b_best_times, double ground_E, double *best_E_array, int seed)
{
	int i=0,j=0, random_row=0, random_col, proposal_accepted=0,proposal_count=0, poor_acceptance_count=0,random_time_index;
	double *psi, *k_times, *j_times,*b_times,*k_times_temp, *j_times_temp, *b_times_temp, acceptance_rate=0, E_old=0, E_new=0,best_E=0, change=0;
	k_times = (double *) malloc(2*NUMBER_OF_BANGS*sizeof(double));
	j_times = (double *) malloc(2*NUMBER_OF_BANGS*sizeof(double));
	b_times = (double *) malloc(2*NUMBER_OF_BANGS*sizeof(double));
	k_times_temp = (double *) malloc(2*NUMBER_OF_BANGS*sizeof(double));
	j_times_temp = (double *) malloc(2*NUMBER_OF_BANGS*sizeof(double));
	b_times_temp = (double *) malloc(2*NUMBER_OF_BANGS*sizeof(double));
	psi = (double *) malloc (2*N*sizeof(double));

	memcpy(psi,psi_start, 2*N*sizeof(double));
	copy_arrays_bang_bang(j_times, k_times, b_times, j_best_times, k_best_times, b_best_times);


	if(CHECK) check_norm(psi, N);

	evolve_mc_bang_bang(table,b,num_electrons,N,lattice, psi, j_times,k_times,b_times, tau);
	best_E = cost(psi, ham_target, N);

	if(tau==TAU_INIT_MCBB) best_E_array[0] = best_E;
	E_old = best_E;

	if(PRINT) printf("Pre-Monte_Carlo Expectation:   %f\n", best_E);







	for (i=0;i<TEMP_DECAY_ITERATIONS_MCBB;i++)
	{
		proposal_accepted = 0;
		proposal_count = 0;

		for (j=0; j<SWEEPS_MCBB*NUMBER_OF_BANGS;j++)
		{

			copy_arrays_bang_bang(j_times_temp, k_times_temp, b_times_temp, j_times, k_times, b_times);
			change = get_change_mcbb(r, tau);


			random_time_index = (int)floor(get_random(0, 2*NUMBER_OF_BANGS, r));
			change_array_bang_bang(j_times,k_times,b_times,change,random_time_index,j, tau);




			memcpy(psi,psi_start, 2*N*sizeof(double));//resetting psi

			evolve_mc_bang_bang(table,b,num_electrons,N,lattice, psi, j_times,k_times,b_times, tau);
			E_new = cost(psi, ham_target, N);
			if (E_new<best_E) best_E=E_new, E_old=E_new,  copy_arrays_bang_bang(j_best_times, k_best_times, b_best_times, j_times, k_times, b_times);
			else if (E_new<=E_old) E_old=E_new;
			else if (get_random(0,1,r)<exp(-(E_new-E_old)/(temp))) E_old=E_new, proposal_accepted++, proposal_count++;
			else copy_arrays_bang_bang( j_times, k_times, b_times,j_times_temp, k_times_temp, b_times_temp), proposal_count++;
			//printf("dif : %f\n\n", E_new-E_old);
		}

		acceptance_rate = (double)proposal_accepted/proposal_count;

		if(PRINT_MC_ITERATIONS) printf("           Best Expectation:   %3.6f  ||  Acceptance Rate: %3.4f (%i/%i)\n",best_E,acceptance_rate,proposal_accepted, proposal_count);

		if(acceptance_rate<0.1) poor_acceptance_count++;
		else poor_acceptance_count = 0;

		if(poor_acceptance_count>TEMP_DECAY_LIMIT_MCBB)
		{
			printf("NO MC PROGRESS FOR %i TEMP_DECAY ITERATIONS, TERMINATING\n", TEMP_DECAY_LIMIT_MCBB);
			break;
		}
		temp=temp*TEMP_EXP_DECAY_MCBB;
	}


	if(PRINTBEST) print_times(j_best_times, k_best_times, b_best_times);
	if(PRINT_TIMES) print_times(j_best_times, k_best_times, b_best_times);

	free(k_times), free(j_times), free(b_times),free(k_times_temp), free(j_times_temp), free(b_times_temp),free(psi);
	return best_E;
}










void adiabatic_method(gsl_rng * r,int *table, unsigned long long int *b, int num_electrons, int N, int lattice[NX][NY], double*ham_target,double* psi_start,double* jkb_initial,double* jkb_target,double g_initial, double f_initial, double g_target, double f_target, double ground_E)
{
	int total_steps, index;
	double *tau_array, *E_array, tau, time_step, *psi;

	tau_array = (double *) malloc(400*sizeof(double));
	E_array = (double *) malloc(400*sizeof(double));
	psi = (double *) malloc (2*N*sizeof(double));

	memcpy(psi,psi_start, 2*N*sizeof(double));

	tau = TAU_INIT_ADIA;
	time_step = TIME_STEP_ADIA;
	total_steps = floor(int(tau/double(time_step)));
	index = 0;

	if(PRINT)
	{
		printf("#######################################################################");
		printf("\n####################### THE ADIABATIC METHOD ##########################\n");
		printf("#######################################################################");
		printf("\n ELECTRONS:        %2i || DIMENSION:    %4i || TAU_MAX:        %4.2f || \n TOTAL_STEPS:    %4i || TIME_STEP:   %4.3f || G_INIT:         %4.3f ||\n F_INIT:        %4.3f || G_TARG:      %4.3f || F_TARG:         %4.3f ||\n",num_electrons, N,double(MAX_TAU_ADIA),total_steps, time_step,g_initial, f_initial, g_target, f_target);
		printf("\nPSI_START: "), print_vector(psi_start, N);
		printf("#######################################################################\n");
	}
	while(tau<MAX_TAU_ADIA)
	{
		memcpy(psi,psi_start, 2*N*sizeof(double));//resetting psi
		evolve_adiabatic(table,b,num_electrons,N,lattice, psi, total_steps, time_step, tau, g_initial, f_initial, g_target, f_target);

		E_array[index] = cost(psi, ham_target, N);
		tau_array[index] = tau;
		printf("       Tau: %5.2f  ||  Expectation-Value:  %7.4f  ||  Target: %f \n", tau, E_array[index],ground_E);

		if((E_array[index] - ground_E) < DIFFERENCE_LIMIT_ADIA)
		{
			printf("\n\nGROUND STATE REACHED. STOPPING THE ADIABATIC EVOLUTION\n\n");
		       	break;
		}
		index++;
		tau=tau+TAU_INIT_ADIA;//TAU_SCALAR_ADIA;
		total_steps = floor(tau/time_step);
	}
	if (ADIABATIC_DATA) export_adiabatic_data(tau_array, E_array, total_steps, tau, ground_E, g_initial, f_initial, g_target, f_target, jkb_initial, jkb_target, index, num_electrons);
	free(tau_array), free(E_array);
}











void evolve_mc(int *table, unsigned long long int *b,int num_electrons,int N, int lattice[NX][NY], double *psi, double *k_array, double *j_array, double *b_array, int *bonds, int total_steps, double time_step)
/*evolve a starting state, psi, by acting on it with exp(ham_dev*-i*time_step). The resulting state is updated as psi and the process repeats total_steps times (tau/total_steps) until the final psi state is produced. The function contains two methods for calculating exp(ham_dev*-i+time_step), one is a diagonalization method, the other a Pade approximation*/
{
	int i,j;
	double *ham_dev,*ham_t_i, *ham_real,*exp_matrix,*D, *Vdag;
	ham_dev = (double *) malloc (2*N*N*sizeof (double));
	exp_matrix = (double *) malloc (2*N*N*sizeof (double));
	ham_t_i = (double *) malloc (2*N*N*sizeof (double));
	ham_real = (double *) malloc(N*N*sizeof(double));
	Vdag = (double*) malloc(N*N*sizeof(double));
	D = (double*) malloc(N*sizeof(double));

	for (i=0; i<total_steps;i++)
	{
		construct_device_hamiltonian(table, b, num_electrons, N, ham_dev, lattice, k_array, j_array, b_array, bonds,i,DEVICE_DIMENSION);

		if(CHECK) check_norm(psi, N);

		if(DIAG)
		{
			for (j=0; j<N*N; j++) ham_real[j]=0.0,Vdag[j]=0.0;
			for (j=0; j<N; j++) D[j]=0.0;
			for (j=0; j<N*N; j++) ham_real[j] = ham_dev[2*j];//converting an all-real-valued complex matrix into just real matrix
			diag_hermitian_real_double(N, ham_real,Vdag, D);
			exp_diaganolized_mat(ham_dev, Vdag, D, N, time_step);//This function exponentiates D to e^(-i*time_step*D)

			if(CHECK) check_unitary(ham_dev, N);
			matrix_vector_mult(ham_dev,psi, N);
		}
		else
		{
			for (j=0; j<N*N*2; j++) exp_matrix[j] = 0.0, ham_t_i[j]=0.0;
			for (j=0; j<N*N; j++)
			{
				ham_t_i[2*j+1] = (ham_dev[2*j]*-time_step); //multiplying by -i*dt for the Pade approximation
				ham_t_i[2*j] = (ham_dev[2*j+1]*time_step);
			}
			exp_general_complex_double(N, ham_t_i, exp_matrix);
			if(CHECK) check_unitary(exp_matrix, N);
			matrix_vector_mult(exp_matrix, psi, N);
		}
		if(CHECK) check_norm(psi, N);
		//print_vector(psi,N);
	}
	free(ham_dev), free(ham_t_i), free(exp_matrix), free(D), free(Vdag),free(ham_real);
}




void evolve_mc_bang_bang(int *table, unsigned long long int *b,int num_electrons,int N,  int lattice[NX][NY], double *psi, double *j_times, double *k_times, double *b_times, double tau)
{
	int i,j, *jkb_index;
	double *ham_dev,*ham_t_i, *ham_real,*exp_matrix,*D, *Vdag, time, time2, time_step, *jkb_vals;

	ham_dev = (double *) malloc (2*N*N*sizeof (double));
	exp_matrix = (double *) malloc (2*N*N*sizeof (double));
	ham_t_i = (double *) malloc (2*N*N*sizeof (double));
	ham_real = (double *) malloc(N*N*sizeof(double));
	Vdag = (double*) malloc(N*N*sizeof(double));
	D = (double*) malloc(N*sizeof(double));
	jkb_index = (int *) malloc(3*sizeof(int));
	jkb_index[0] = 0;
	jkb_index[1] = 0;
       	jkb_index[2] = 0;
	jkb_vals = (double *) malloc(3*sizeof(double));
	time = 0;

	if(PRINT_TIMES) print_times(j_times, k_times, b_times);

	while(time < tau)
	{

		//printf("jkb_index = %i, %i, %i\n", jkb_index[0],jkb_index[1],jkb_index[2]);
		time2 = find_next_time(time, tau,j_times, k_times, b_times, jkb_index, jkb_vals);
		time_step = time2-time;
		//printf("time: %f, timestep = %f\n", time, time_step);
		time = time2;

		//printf("jkb : %f, %f, %f\n", jkb[0],jkb[1],jkb[2]);
		construct_device_hamiltonian_uniform(table, b, num_electrons, N, ham_dev, lattice, jkb_vals,DEVICE_DIMENSION);
		if(CHECK) check_norm(psi, N);

		if(DIAG)
		{
			for (j=0; j<N*N; j++) ham_real[j]=0.0,Vdag[j]=0.0;
			for (j=0; j<N; j++) D[j]=0.0;
			for (j=0; j<N*N; j++) ham_real[j] = ham_dev[2*j];//converting an all-real-valued complex matrix into just real matrix
			diag_hermitian_real_double(N, ham_real,Vdag, D);
			exp_diaganolized_mat(ham_dev, Vdag, D, N, time_step);//This function exponentiates D to e^(-i*time_step*D)

			if(CHECK) check_unitary(ham_dev, N);
			matrix_vector_mult(ham_dev,psi, N);
		}
		else
		{
			for (j=0; j<N*N*2; j++) exp_matrix[j] = 0.0, ham_t_i[j]=0.0;
			for (j=0; j<N*N; j++)
			{
				ham_t_i[2*j+1] = (ham_dev[2*j]*-time_step); //multiplying by -i*dt for the Pade approximation
				ham_t_i[2*j] = (ham_dev[2*j+1]*time_step);
			}
			exp_general_complex_double(N, ham_t_i, exp_matrix);
			if(CHECK) check_unitary(exp_matrix, N);
			matrix_vector_mult(exp_matrix, psi, N);
		}
		if(CHECK) check_norm(psi, N);
	}
	free(ham_dev), free(ham_t_i), free(exp_matrix), free(D), free(Vdag),free(ham_real), free(jkb_vals), free(jkb_index);
}



void print_times(double* j_times, double* k_times, double* b_times)
{
	int i;
	printf("\nj_times:");

	for (i=0; i<2*NUMBER_OF_BANGS; i++) printf(" %5.4f |", j_times[i]);

	printf("\nk_times:");
	for (i=0; i<2*NUMBER_OF_BANGS; i++) printf(" %5.4f |", k_times[i]);

	printf("\nb_times:");
	for (i=0; i<2*NUMBER_OF_BANGS; i++) printf(" %5.4f |", b_times[i]);
	printf("\n");
}





void evolve_adiabatic(int *table, unsigned long long int *b,int num_electrons,int N, int lattice[NX][NY], double *psi, int total_steps, double time_step,double tau, double g_i, double f_i, double g_t, double f_t)
/*evolve a starting state, psi, by acting on it with exp(ham_dev*-i*time_step). The resulting state is updated as psi and the process repeats total_steps times (tau/total_steps) until the final psi state is produced. The function contains two methods for calculating exp(ham_dev*-i+time_step), one is a diagonalization method, the other a Pade approximation*/
{
	int i,j;
	double *ham,*ham_t_i, *ham_real,*exp_matrix,*D, *Vdag, g, f, *jkb;
	ham = (double *) malloc (2*N*N*sizeof (double));
	exp_matrix = (double *) malloc (2*N*N*sizeof (double));
	ham_t_i = (double *) malloc (2*N*N*sizeof (double));
	ham_real = (double *) malloc(N*N*sizeof(double));
	Vdag = (double*) malloc(N*N*sizeof(double));
	D = (double*) malloc(N*sizeof(double));
	jkb = (double *) malloc(3*sizeof(double));

	for (i=1; i<total_steps+1;i++)
	{

		g = g_i + (g_t-g_i)*((i*1.0)/total_steps);
		f = f_i + (f_t-f_i)*(i/total_steps);
		jkb[0] = g;
		jkb[1] = (1-g)*f;
		jkb[2] = (1-f)*(1-g);

		construct_device_hamiltonian_uniform(table, b, num_electrons, N, ham, lattice, jkb, DEVICE_DIMENSION);

		if(CHECK) check_norm(psi, N);

		if(DIAG)
		{
			for (j=0; j<N*N; j++) ham_real[j]=0.0,Vdag[j]=0.0;
			for (j=0; j<N; j++) D[j]=0.0;
			for (j=0; j<N*N; j++) ham_real[j] = ham[2*j];//converting an all-real-valued complex matrix into just real matrix
			diag_hermitian_real_double(N, ham_real,Vdag, D);
			exp_diaganolized_mat(ham, Vdag, D, N, time_step);//This function exponentiates D to e^(-i*time_step*D)

			if(CHECK) check_unitary(ham, N);
			matrix_vector_mult(ham,psi, N);
		}
		else
		{
			for (j=0; j<N*N*2; j++) exp_matrix[j] = 0.0, ham_t_i[j]=0.0;
			for (j=0; j<N*N; j++)
			{
				ham_t_i[2*j+1] = (ham[2*j]*-time_step); //multiplying by -i*dt for the Pade approximation
				ham_t_i[2*j] = (ham[2*j+1]*time_step);
			}
			exp_general_complex_double(N, ham_t_i, exp_matrix);
			if(CHECK) check_unitary(exp_matrix, N);
			matrix_vector_mult(exp_matrix, psi, N);
		}
		if(CHECK) check_norm(psi, N);
	}
	free(ham), free(ham_t_i), free(exp_matrix), free(D), free(Vdag),free(ham_real);
}






double cost(double *psi, double *ham_target, int N)
/*Computing the expectation value between ham_target and psi, <psi|ham_target|psi>*/
{

	int i=0,INCX = 1,INCY = 1,j;
	double *psi_conj, *psi_temp, result=0, resulti=0, norm;
	psi_conj = (double*) malloc (N*2*sizeof(double));
	psi_temp = (double*) malloc (N*2*sizeof(double));
	memcpy(psi_conj, psi, 2*N*sizeof(double));
	memcpy(psi_temp, psi, 2*N*sizeof(double));


	if(CHECK) check_norm(psi_temp, N);
	if(CHECK) check_weights(psi_temp, ham_target, N);

	matrix_vector_mult(ham_target, psi_temp, N);//H*psi, the operator acting on the ket, storing in psi
	result = zdotc_(&N, psi_conj, &INCX, psi_temp, &INCY);//psi* * psi, the dot between the complex conj and the result of H*psi

	free(psi_conj), free(psi_temp);
	return result;
}




double calc_initial_temp(int *table, unsigned long long int *b, int num_electrons, int N, int lattice[NX][NY], double*ham_target, int total_steps, double time_step,double tau, double *psi_start, gsl_rng *r, int seed, double *jkb_initial, double *jkb_target)
/*Finding an good initial temp that allows an average increase acceptance probability of about ACCEPTANCE_PROB_MC (typically 0.8). Choosing an initial temperature that allows the chance of accepting j,b, and k_array values which increase the expectation value to be ~80%, https://www.phy.ornl.gov/csep/mo/node32.html*/
{
	if(PRINT) printf("\n\n\n\n...Calculating initial temperature based on %i random starting states...\n", RANDOM_STATES_MC);
	int *bonds, i=0,j=0, random_row=0,random_col=0, start_state=0, count=0;
	double *psi,*psi_random, *k_array, *j_array,*b_array, E_old=0, E_new=0,sum=0, change, initial_temp=0;
	bonds = (int*) malloc(NUM_SITES*NUM_SITES*sizeof(int));
	k_array = (double *) malloc(2*NUM_SITES*total_steps*sizeof(double));
	j_array = (double *) malloc(2*NUM_SITES*total_steps*sizeof(double));
	b_array = (double *) malloc(NUM_SITES*total_steps*sizeof(double));
	psi_random = (double *) malloc (2*N*sizeof(double));
	psi = (double *) malloc (2*N*sizeof(double));

	assign_bonds(bonds, lattice);
	for (i=0; i<N*2;i++) psi[i]=0.0;
	for (j=0;j<RANDOM_STATES_MC;j++)
	{
		for (i=0; i<N*2;i++) psi_random[i] =0.0;
		start_state = floor(get_random(0,N,r));

		psi_random[start_state*2] = 1;
		memcpy(psi,psi_random, 2*N*sizeof(double));

		init_arrays(k_array, j_array, b_array,r, total_steps);

		evolve_mc(table,b,num_electrons,N,lattice, psi, k_array,j_array,b_array, bonds, total_steps, time_step);
		E_old = cost(psi, ham_target, N);

		for (i=0; i<SWEEPS_MC*total_steps;i++)
		{
			change = get_change_mc(r, tau);
			random_row = floor(get_random(0,total_steps,r));
			random_col = floor(get_random(0,NUM_SITES*2,r));

			change_array(k_array,j_array,b_array,random_row,random_col,change,i, total_steps);
			memcpy(psi,psi_random, 2*N*sizeof(double));//resetting psi
			evolve_mc(table,b,num_electrons,N,lattice, psi, k_array,j_array,b_array, bonds, total_steps, time_step);
			E_new = cost(psi, ham_target, N);

			if (E_new>E_old) sum += (E_new-E_old), count++;
			E_old=E_new;
		}
	}
	free(k_array), free(j_array), free(b_array), free(psi), free(psi_random), free(bonds);
	initial_temp = -(sum/(count*log(ACCEPTANCE_PROB_MC)));

	if(PRINT)
	{
		printf("#######################################################################");
		printf("\n################# USING THE MONTE-CARLO METHOD ########################");
		printf("\n#######################################################################");
		printf("\n ELECTRONS:         %2i || DIMENSION:      %4i || SEED:         %4i ||\n TAU:           %6.4f || TOTAL_STEPS:    %4i || TIME_STEP:   %4.3f ||\n INIT_TEMP:     %5.4f || TOTAL_SWEEPS:   %4i || TEMP_DECAYS:  %4i ||\n J_INIT:         %4.3f || K_INIT:        %4.3f || B_INIT:      %4.3f ||\n J_TARGET:       %4.3f || K_TARGET:      %4.3f || B_TARGET:    %4.3f ||\n",num_electrons, N,seed,time_step*total_steps,total_steps, time_step,initial_temp, total_steps*SWEEPS_MC, TEMP_DECAY_ITERATIONS_MC, jkb_initial[0], jkb_initial[1], jkb_initial[2], jkb_target[0], jkb_target[1], jkb_target[2]);
		printf("\nPSI_START: "), print_vector(psi_start, N);
		printf("#######################################################################\n");
	}

	return initial_temp;
}



double calc_initial_temp_bang_bang(int *table, unsigned long long int *b, int num_electrons, int N, int lattice[NX][NY], double*ham_target, double tau, double *psi_start, gsl_rng *r, int seed, double *jkb_initial, double *jkb_target)
/*Finding an good initial temp that allows an average increase acceptance probability of about ACCEPTANCE_PROB_MC (typically 0.8). Choosing an initial temperature that allows the chance of accepting j,b, and k_array values which increase the expectation value to be ~80%, https://www.phy.ornl.gov/csep/mo/node32.html*/
{

	if(PRINT) printf("\n\n\n\n...Calculating initial temperature based on %i random starting states...\n", RANDOM_STATES_MCBB);

	int i=0,j=0, random_row=0,random_col=0, start_state=0, count=0;
	double *psi,*psi_random, *k_times;
	double *j_times;
	double *b_times;
        double E_old=0;
	double E_new=0,sum=0, change=0, initial_temp=0,random_time_index;


	k_times = (double *) malloc(2*NUMBER_OF_BANGS*sizeof(double));
	j_times = (double *) malloc(2*NUMBER_OF_BANGS*sizeof(double));
	b_times = (double *) malloc(2*NUMBER_OF_BANGS*sizeof(double));

	psi_random = (double *) malloc (2*N*sizeof(double));
	psi = (double *) malloc (2*N*sizeof(double));

	for (i=0; i<N*2;i++) psi[i]=0.0;
	for (j=0;j<RANDOM_STATES_MCBB;j++)
	{
		for (i=0; i<N*2;i++) psi_random[i] =0.0;
		start_state = floor(get_random(0,N,r));

		psi_random[start_state*2] = 1;
		memcpy(psi,psi_random, 2*N*sizeof(double));

		init_arrays_bang(j_times, k_times, b_times, r);

		evolve_mc_bang_bang(table,b,num_electrons,N,lattice, psi, j_times,k_times,b_times, tau);
		E_old = cost(psi, ham_target, N);

		for (i=0; i<SWEEPS_MCBB*NUMBER_OF_BANGS;i++)
		{
			change = get_change_mcbb(r, tau);
			random_time_index = (int)floor(get_random(0, 2*NUMBER_OF_BANGS, r));
			change_array_bang_bang(j_times,k_times,b_times,change,random_time_index,j, tau);
			memcpy(psi,psi_start, 2*N*sizeof(double));//resetting psi
			evolve_mc_bang_bang(table,b,num_electrons,N,lattice, psi, j_times,k_times,b_times, tau);

			E_new = cost(psi, ham_target, N);
			if (E_new>=E_old) sum += (E_new-E_old), count++;
			E_old=E_new;
		}
	}
	free(k_times), free(j_times), free(b_times), free(psi), free(psi_random);
	initial_temp = -(sum/(count*log(ACCEPTANCE_PROB_MC)));

	if(PRINT)
	{
		printf("#######################################################################");
		printf("\n################# USING THE MC BANG-BANG METHOD #######################");
		printf("\n#######################################################################");
		printf("\n ELECTRONS:         %2i || DIMENSION:      %4i || SEED:         %4i ||\n TAU:           %6.4f || NUMBER OF BANGS:%4i ||                    ||\n INIT_TEMP:     %5.4f || TOTAL_SWEEPS:   %4i || TEMP_DECAYS:  %4i ||\n J_INIT:         %4.3f || K_INIT:        %4.3f || B_INIT:      %4.3f ||\n J_TARGET:       %4.3f || K_TARGET:      %4.3f || B_TARGET:    %4.3f ||\n",num_electrons, N,seed,tau,NUMBER_OF_BANGS,initial_temp, SWEEPS_MCBB*NUMBER_OF_BANGS, TEMP_DECAY_ITERATIONS_MCBB, jkb_initial[0], jkb_initial[1], jkb_initial[2], jkb_target[0], jkb_target[1], jkb_target[2]);
		printf("\nPSI_START: "), print_vector(psi_start, N);
		printf("#######################################################################\n");
	}

	return initial_temp;
}



double get_change_mc(gsl_rng * r, double tau)
{
	int sign = pow(-1,(int)floor(get_random(0,10,r)));
	double max_change = MAX_CHANGE_MC_INIT * (TAU_INIT_MC/tau);
       	double abs_change = get_random(0,max_change, r);
	return sign*abs_change;
}

double get_change_mcbb(gsl_rng * r, double tau)
{
	int sign = pow(-1,(int)floor(get_random(0,10,r)));
	double max_change = tau*CHANGE_FRACTION_MCBB;
       	double abs_change = get_random(0,max_change, r);
	return sign*abs_change;
}



void construct_device_hamiltonian(int *table, unsigned long long int *b,int num_electrons,int N, double *ham_dev, int lattice[NX][NY],double *k_array, double *j_array, double *b_array, int *bonds,  int index, int D)
/*constructiing the hamiltonina matrix for the device amiltonian, using the index-ith row of each j, b, and k array*/
{
	int i=0,ii=0,j=0,x=0,y=0,state=0,site=0,sign=0,bond=0,neighbor_count=0,*neighbors;
	unsigned long long int *v,comparison=0;
	v = (unsigned long long int*) malloc(sizeof(unsigned long long int));
	neighbors = (int*) malloc(4*sizeof(int));
	for (i=0; i<N*N*2; i++) ham_dev[i] = 0.0;

	for (i=0;i<N;i++)
	{

		for (j=0;j<num_electrons;j++)//The J term calculation
		{

			site=table[i*num_electrons+j];
			neighbor_count = get_neighbors(site, neighbors, lattice, D);
			for (ii=0; ii<neighbor_count; ii++)
			{
				if (((1ULL<<(neighbors[ii]-1))&b[i])==0)//making sure neighbor is not occupied, otherwise nothing happens
				{

					hop(b[i], v,site, neighbors[ii]);

					state=find(N,v, b);

					bond = bonds[(site-1)*NUM_SITES+neighbors[ii]-1]-1;
					ham_dev[(N*i+state-1)*2] -= j_array[index*NUM_SITES*2+bond];
				}
			}
		}

		for (j=1;j<(NUM_SITES);j++)//The K term calculation
		{
			site=j;
			neighbor_count = get_neighbors(site, neighbors, lattice, D);
			for (ii=0; ii<neighbor_count;ii++)
			{
				if (neighbors[ii] > site)
				{
					sign = -1;
					comparison = (1ULL<<(neighbors[ii]-1))+(1ULL<<(site-1));
					if((comparison&b[i])==comparison || (comparison&b[i])==0) sign = 1;
					bond = bonds[(site-1)*NUM_SITES+neighbors[ii]-1]-1;
					ham_dev[((N*i)+i)*2] += k_array[index*NUM_SITES*2+bond]*sign;
				}
			}
		}

		for (j=0; j<NUM_SITES;j++)//The B term calculation
		{
			sign = -1;
			if(((1ULL<<j)&b[i])>0) sign=1;
			ham_dev[((N*i)+i)*2] += b_array[index*NUM_SITES+j]*sign;
		}
	}
	b = (unsigned long long int*) malloc(N*sizeof(double));
	free(neighbors);// free(v);
	if(CHECK) check_hermicity(ham_dev, N);
}



void construct_device_hamiltonian_uniform(int *table, unsigned long long int *b,int num_electrons,int N, double *ham_dev, int lattice[NX][NY],double *jkb, int D)
/*constructiing the hamiltonina matrix for the device amiltonian, using the index-ith row of each j, b, and k array*/
{
	int i=0,ii=0,j=0,x=0,y=0,state=0,site=0,sign=0,neighbor_count=0,*neighbors;
	unsigned long long int *v,comparison=0;
	v = (unsigned long long int*) malloc(sizeof(unsigned long long int));
	neighbors = (int*) malloc(4*sizeof(int));
	for (i=0; i<N*N*2; i++) ham_dev[i] = 0.0;

	for (i=0;i<N;i++)
	{
		for (j=0;j<num_electrons;j++)//The J term calculation
		{

			site=table[i*num_electrons+j];
			neighbor_count = get_neighbors(site, neighbors, lattice, D);


			for (ii=0; ii<neighbor_count; ii++)
			{

				if (((1ULL<<(neighbors[ii]-1))&b[i])==0)//making sure neighbor is not occupied, otherwise nothing happens
				{
					hop(b[i], v,site, neighbors[ii]);
					state=find(N,v, b);
					ham_dev[(N*i+state-1)*2] -= jkb[0];
				}
			}
		}

		for (j=1;j<(NUM_SITES);j++)//The K term calculation
		{
			site=j;
			neighbor_count = get_neighbors(site, neighbors, lattice, D);
			for (ii=0; ii<neighbor_count;ii++)
			{
				if (neighbors[ii] > site)
				{
					sign = -1;
					comparison = (1ULL<<(neighbors[ii]-1))+(1ULL<<(site-1));
					if((comparison&b[i])==comparison || (comparison&b[i])==0) sign = 1;
					ham_dev[((N*i)+i)*2] += jkb[1]*sign;
				}
			}
		}

		for (j=0; j<NUM_SITES;j++)//The B term calculation
		{
			sign = -1;
			if(((1ULL<<j)&b[i])>0) sign=1;
			ham_dev[((N*i)+i)*2] += jkb[2]*sign;
		}
	}
	free(neighbors), free(v);
	if(CHECK) check_hermicity(ham_dev, N);
}


void construct_model_hamiltonian(int *table, unsigned long long int *b,int num_electrons,int N, double *ham_mod, int lattice[NX][NY], int D)
/*Constructing the hamiltonian matrix for the model hamiltonian*/
{
	int i,ii,j,x,y,state,site,neighbor,sign,neighbor_count, *neighbors;
	unsigned long long int *v,comparison;
	v = (unsigned long long int*) malloc(sizeof(unsigned long long int));
	neighbors = (int*) malloc(sizeof(int)*4);
	for (i=0; i<N*N*2; i++) ham_mod[i] = 0;

	for (i=0;i<N;i++)
	{
		for (j=0;j<num_electrons;j++)//The T term calculation
		{
			site=table[i*num_electrons+j];

			neighbor_count = get_neighbors(site, neighbors, lattice, D);
			for (ii=0; ii<neighbor_count; ii++)
			{
				if (((1ULL<<(neighbors[ii]-1))&b[i])==0)//checking if the neighbor is occupied
				{
					sign = hop(b[i], v,site, neighbors[ii]);
					if (sign==0) sign=1;
					else sign=-1;
					state=find(N,v, b);
					ham_mod[(N*i+state-1)*2] -= (T*sign);
				}
			}
		}
		for (j=1;j<(NUM_SITES);j++)//The V term calculation
		{
			site=j;
			neighbor_count = get_neighbors(site, neighbors, lattice, D);
			for (ii=0; ii<neighbor_count;ii++)
			{
				if (neighbors[ii] > site)
				{
					sign = -1;
					comparison = (1ULL<<(neighbors[ii]-1))+(1ULL<<(site-1));
					if((comparison&b[i])==comparison || (comparison&b[i])==0) sign = 1;
					ham_mod[((N*i)+i)*2] += sign*V;
				}
			}
		}
	}
	free(neighbors), free(v);
	if(CHECK) check_hermicity(ham_mod, N);
}


double find_next_time(double time, double tau,double* j_times, double* k_times,double*  b_times,int*  jkb_index, double*jkb)
{

	bool j, k, b;
	double newtime=tau;
	j = true, k = true, b =true;

	jkb[0] = fmod(jkb_index[0],2), jkb[1] = fmod(jkb_index[1],2), jkb[2] = fmod(jkb_index[2],2);
	if(jkb_index[0] == NUMBER_OF_BANGS*2) j =false;
	if(jkb_index[1] == NUMBER_OF_BANGS*2) k =false;
	if(jkb_index[2] == NUMBER_OF_BANGS*2) b =false;

	if(!(j || k ||b)) return tau;
	//if(j && j_times[jkb_index[0]] == 0.0) jkb_index[0] += 1;
	//if(k && k_times[jkb_index[1]] == 0.0) jkb_index[1] += 1;
	//if(b && b_times[jkb_index[2]] == 0.0) jkb_index[2] += 1;


	if(j) newtime = j_times[jkb_index[0]];
	if(k  && (k_times[jkb_index[1]] < newtime)) newtime = k_times[jkb_index[1]];
	if(b  && (b_times[jkb_index[2]] < newtime)) newtime = b_times[jkb_index[2]];
	if(newtime == j_times[jkb_index[0]]) jkb_index[0] += 1;
	if(newtime == k_times[jkb_index[1]]) jkb_index[1] += 1;
	if(newtime == b_times[jkb_index[2]]) jkb_index[2] += 1;

	if(jkb_index[0] > 2*NUMBER_OF_BANGS ||  jkb_index[0] > 2*NUMBER_OF_BANGS || jkb_index[0] > 2*NUMBER_OF_BANGS) printf("ERROR: FIND_NEXT_TIME INDEX OUT OF BOUNDS"), exit(0);

/*


	bool j, k, b;
	double newtime=0;
	j = true, k = true, b =true;
       	if(jkb_index[0] == 2*NUMBER_OF_BANGS-1) newtime = j_times[jkb_index[0] + 1], j =false;
	if (jkb_index[1] == 2*NUMBER_OF_BANGS-1) newtime = k_times[jkb_index[1] + 1], k=false;
	if (jkb_index[2] == 2*NUMBER_OF_BANGS-1) newtime = b_times[jkb_index[2] + 1], b=false;
	if(!(j || k || b)) newtime = tau;//, printf("WHAT\n\n");
 *
	if(j)
	{
		newtime = j_times[jkb_index[0]+1];
	}
	if(k && k_times[jkb_index[1] + 1] <= newtime)
	{
		newtime = k_times[jkb_index[1]+1];
	}
	if(b && b_times[jkb_index[2] + 1] <= newtime)
	{
		newtime = b_times[jkb_index[2]+1];
	}
	//printf("newtime: %f\n\n", newtime);
	if (newtime == j_times[jkb_index[0]+1])
	{
		jkb_index[0] +=1;
		if (jkb[0] == 1.0) jkb[0] = 0;
		else jkb[0] = 1;
	  //     	printf("j:%f\n\n",jkb[0]);
	}
	if (newtime == k_times[jkb_index[1]+1])
	{
	       	jkb_index[1] +=1;
		if (jkb[1] == 1.0) jkb[1] = 0;
		else jkb[1] = 1.0;
	    //   	printf("k:%f\n\n",jkb[1]);
	}
	if (newtime == b_times[jkb_index[2]+1])
	{
		jkb_index[2] += 1;
		if (jkb[2] == 1.0) jkb[2] = 0;
		else jkb[2] = 1;
	//	printf("b:%f\n\n",jkb[2]);
	}*/
//	printf("jkb: %f %f %f\n", jkb[0], jkb[1], jkb[2]);
//	printf("NEWTIME: %f | j_ind: %i | k_ind: %i | b_ind: %i\n", newtime, jkb_index[0], jkb_index[1], jkb_index[2]);
	return newtime;
}



void diag_hermitian_real_double(int N,  double *A, double *Vdag,double *D)
/*diagonalizing an real square matrix A. Stores the eigenvectors in Vdag and the eigenvalues in D*/
{
	char JOBZ='V',UPLO='U';
	int LDA=N,LWORK=-1, INFO;
	double *WORK;
	WORK=(double*) malloc(sizeof(double));

	memcpy(Vdag,A,N*N*sizeof(double));

	dsyev_(&JOBZ, &UPLO, &N, Vdag, &LDA, D, WORK, &LWORK, &INFO );
	if (INFO !=0) printf("\n\n\nDIAGONALIZATION ERROR, INFO = %i\n\n\n", INFO);
       	LWORK=WORK[0];
       	free(WORK);
       	WORK=(double*) malloc(LWORK*sizeof(double));

       	dsyev_(&JOBZ, &UPLO, &N, Vdag, &LDA, D, WORK, &LWORK, &INFO );
       	if (INFO !=0) printf("\n\n\nDIAGONALIZATION ERROR, INFO = %i\n\n\n", INFO);
       	free(WORK);
}




void exp_diaganolized_mat(double *ham, double *Vdag, double* D, int N, double time_step)
/*calculating the exponential of a diagonalized decomposition where the matrix A = Vdag*D*Vdag_inv. calculating Vdag*exp(D)*Vdag_inv =exp(A), storing the result in ham*/
{
	char TRANSA = 'N', TRANSB = 'N';
	int i, *IPIV, LWORK=N*N, INFO, LDA=N, LDB=N, LDC=N;
	double *exp_D, *temp_mat, *Vdag_z, *Vdag_z_inv, *WORK,ALPHA[2], BETA[2];
	ALPHA[0]=1.0, ALPHA[1]=0.0;
	BETA[0]=0.0, BETA[1]=0.0;

	IPIV = (int*) malloc(N*sizeof(int));
	WORK = (double*) malloc(LWORK*sizeof(double));
	exp_D = (double *) malloc (2*N*N*sizeof (double));
	temp_mat = (double *) malloc (2*N*N*sizeof (double));
	Vdag_z = (double *) malloc (2*N*N*sizeof (double));
	Vdag_z_inv = (double *) malloc (2*N*N*sizeof (double));

	for (i =0; i<N*N*2; i++) exp_D[i] = 0, temp_mat[i] = 0,Vdag_z[i]=0, Vdag_z_inv[i]=0, ham[i]=0;
	for (i =0; i<N*N; i++) Vdag_z[2*i] = Vdag[i];
	for (i =0; i<N; i++)
	{
		exp_D[2*(i*N+i)] = cos(-time_step*D[i]);
	       	exp_D[2*(i*N+i)+1] = sin(-time_step*D[i]);
	}

	dgetrf_(&N,&N,Vdag,&N,IPIV,&INFO);
	dgetri_(&N,Vdag,&N,IPIV,WORK,&LWORK,&INFO);//inverting vdag

	for (i =0; i<N*N; i++) Vdag_z_inv[2*i] = Vdag[i];
	zgemm_(&TRANSA, &TRANSB, &N, &N, &N, ALPHA, Vdag_z, &LDA, exp_D, &LDB, BETA, temp_mat, &LDC); //matrix mult
	zgemm_(&TRANSA, &TRANSB, &N, &N, &N, ALPHA, temp_mat, &LDA, Vdag_z_inv, &LDB, BETA, ham, &LDC); //matrix mult
	free(exp_D),free(temp_mat),free(Vdag_z), free(Vdag_z_inv), free(IPIV), free(WORK);
}




void exp_general_complex_double(int N, double *A, double *B)
/*Using the Pade Approximation to calulate exp(A), where A is a complex matrix storing real values at even indices and their imaginary parts at index +1. Storing the result in B*/
{
	int M=N,K=N,ii,jj,kk,s,p,q,INFO,LDA=N,LDB=N,LDC=N,NRHS=N, *IPIV;
	double *row_norm,*X,*Y,*Z,*E,*D,norm,c, ALPHA[2], BETA[2];
	char TRANSA='N',TRANSB='N',TRANS='N';
	ALPHA[0]=1.0,ALPHA[1]=0.0;
	BETA[0]=0.0,BETA[1]=0.0;

	row_norm=(double*) malloc(N*sizeof(double));
	X=(double*) malloc(2*N*N*sizeof(double));
	Y=(double*) malloc(2*N*N*sizeof(double));
	Z=(double*) malloc(2*N*N*sizeof(double));
	E=(double*) malloc(2*N*N*sizeof(double));
	D=(double*) malloc(2*N*N*sizeof(double));
	IPIV=(int*) malloc(N*sizeof(int));
	memcpy(Z,A,2*N*N*sizeof(double));

	for(ii=0;ii<N;ii++)
	{
		row_norm[ii]=0.0;
		for (jj=0;jj<N;jj++) row_norm[ii]=row_norm[ii]+cabs(A[2*(jj+N*ii)]+I*A[2*(jj+N*ii)+1]);

	}
	norm=row_norm[0];
	for(ii=1;ii<N;ii++) if (row_norm[ii]>norm) norm=row_norm[ii];
	s=(int) floor(log2(norm)+2);
	if (s<0) s=0;

	for(ii=0;ii<2*N*N;ii++) Z[ii]=Z[ii]/pow(2.0,s);

	memcpy(X,Z,2*N*N*sizeof(double));
	c=0.5,p=1,q=6;

	for(ii=0;ii<2*N*N;ii++) E[ii]=c*Z[ii];
	for(ii=0;ii<2*N*N;ii++) D[ii]=-c*Z[ii];
	for(ii=0;ii<N;ii++)E[2*(ii+N*ii)]=E[2*(ii+N*ii)]+1.0;
	for(ii=0;ii<N;ii++)D[2*(ii+N*ii)]=D[2*(ii+N*ii)]+1.0;

	for (kk=2;kk<=q;kk++)
	{
		c = c * (q-kk+1) / (kk*(2*q-kk+1));
		zgemm_ (&TRANSA, &TRANSB, &M, &N, &K, ALPHA, Z, &LDA, X, &LDB, BETA, Y, &LDC); //matrix mult
		memcpy(X,Y,2*N*N*sizeof(double));
		for(ii=0;ii<2*N*N;ii++) Y[ii]=c*X[ii];
		for(ii=0;ii<2*N*N;ii++) E[ii]=E[ii]+Y[ii];
		if (p==1)for(ii=0;ii<2*N*N;ii++) D[ii]=D[ii]+Y[ii];
		else for(ii=0;ii<2*N*N;ii++) D[ii]=D[ii]-Y[ii];
		p=abs(1-p);
	}

	zgetrf_ (&M, &N, D, &LDA, IPIV, &INFO);//getting a factorization of D
	if (INFO !=0) printf("\n\n\nERROR, INFO = %i\n\n\n", INFO);
	zgetrs_( &TRANS, &N, &NRHS,D,    &LDA, IPIV,E, &LDB, &INFO);//solving a system of equations
	if (INFO !=0) printf("\n\n\nERROR, INFO = %i\n\n\n", INFO);

	for (kk=1;kk<=s;kk++)
	{
		memcpy(X,E,2*N*N*sizeof(double));
		zgemm_ (&TRANSA, &TRANSB, &M, &N, &K, ALPHA, E, &LDA, X, &LDB, BETA, Y, &LDC);//matrixmultiplication
		if (INFO !=0) printf("\n\n\nERROR, INFO = %i\n\n\n", INFO);
		memcpy(E,Y,2*N*N*sizeof(double));
	}
	memcpy(B,E,2*N*N*sizeof(double));

	free(row_norm),free(X),free(Y),free(Z),free(E),free(D),free(IPIV);
}




void matrix_vector_mult(double *matrix, double *psi, int N)
/*Matrix vector multiplication, psi=matrix*psi*/
{
	char TRANS = 'N';
	int i,INCX = 1,INCY = 1,LDA = N,M=N;
	double *result, ALPHA[2], BETA[2];
	ALPHA[0]=1.0,ALPHA[1]=0.0;
	BETA[0]=0.0,BETA[1]=0.0;
	result = (double *) malloc (N*2*sizeof(double));

	zgemv_(&TRANS, &M, &N,ALPHA,matrix, &LDA, psi, &INCX, BETA, result, &INCY);
	memcpy(psi,result, 2*N*sizeof(double));

	free(result);
}




int hop(unsigned long long int b, unsigned long long int *v,int n, int j)
/*given a state b generates the state v obtained by hopping from site b[n] to site j (which should not be in b), and outputs the fermionic sign*/
{
	unsigned long long int i,x,y;
	int z_count = 0;

//	printf("site %i\n", n);
	//printf("neihbor %i\n", j);
//	printf("b %llu\n",b);


	x = (1ULL << (n-1)) + (1ULL << (j-1));
	for (i=n;i<j-1;i++)  if((1ULL<<i) & (b)) z_count++;
	y = (x ^ b);
	memcpy(v, &y, sizeof(unsigned long long int));
	return z_count%2;
}




int find(int N,unsigned long long int *v,unsigned long long int *b)
/*find the position of a given combination v (vector of length k) in the table of all combinations tab*/
{
	int first, last, mid;
	first=0;
	last=N-1;


	while (first <= last)
	{
		mid = (int) ((first + last) / 2.0);
		if (*v > b[mid]) first = mid + 1;
		else if (*v < b[mid]) last = mid - 1;
		else return mid+1;
	}
}





void check_commutator(int N, double* A, double* B)
{
	char TRANSA = 'N', TRANSB = 'N';
	int i, *IPIV, LWORK=N*N, INFO, LDA=N, LDB=N, LDC=N;
	double *C, *AB, *BA ,ALPHA[2], BETA[2];
	bool commute = true;
	ALPHA[0]=1.0, ALPHA[1]=0.0;
	BETA[0]=0.0, BETA[1]=0.0;


	C = (double*) malloc(2*N*N*sizeof(double));
	BA = (double*) malloc(2*N*N*sizeof(double));
	AB = (double*) malloc(2*N*N*sizeof(double));
	for (i =0; i<N*N*2; i++) C[i] = 0, AB[i] = 0, BA[i]=0;

	zgemm_(&TRANSA, &TRANSB, &N, &N, &N, ALPHA, A, &LDA, B, &LDB, BETA, AB, &LDC); //matrix mult
	zgemm_(&TRANSA, &TRANSB, &N, &N, &N, ALPHA, B, &LDA, A, &LDB, BETA, BA, &LDC); //matrix mult
	for (i =0; i<N*N*2; i++) C[i] = AB[i] - BA[i];
	for (i =0; i<N*N*2; i++) if(C[i] < -0.001 || 0.001 < C[i]) commute = false;
	if(commute) printf("\n\n\nWARNING: THE TARGET AND INITIAL HAMILTONIAN COMMUTE\n\n\n");
	if(PRINT_COMMUTATOR) print_hamiltonian(C, N);

//	print_hamiltonian(B, N);
	free(C), free(BA), free(AB);
}







void get_ground_state(int N, double *ham, double *ground_state)
/*Get the ground state of a hamiltonian matrix by finding the smallest eigen_value, the getting the eigenvector associated with that eigenvalue. Copying the eigenvector the ground_state*/
{
	int i, min_i,j;
	double *D, *Vdag, *ham_real, min_E=1000;
	Vdag = (double*) malloc(N*N*sizeof(double));
	D = (double*) malloc(N*sizeof(double));
	ham_real = (double*) malloc(N*N*sizeof(double));
	for(i=0;i<N*N;i++) ham_real[i] = ham[2*i];

	diag_hermitian_real_double(N, ham_real,Vdag, D);


	for(i=0; i<2*N; i++) ground_state[i] = 0.0;
	for(i=0; i<N; i++)
	{	if(D[i]<min_E)
		{
			min_E = D[i];
			min_i = i;
		}
	}
	for(i=0; i<N; i++) ground_state[i*2] = Vdag[(min_i)*N+i];

	free(Vdag), free(D), free(ham_real);
}




double get_ground_E(int N, double *ham)
/*Finding the ground energy, the smallest eigenvalue of the matrix*/
{
	int i,j;
	double *D, *Vdag, *ham_real,min_E=1000;
	Vdag = (double*) malloc(N*N*sizeof(double));
	D = (double*) malloc(N*sizeof(double));
	ham_real = (double*) malloc(N*N*sizeof(double));
	for(i=0;i<N*N;i++) ham_real[i] = ham[2*i];

	diag_hermitian_real_double(N, ham_real,Vdag, D);

	for(i=0; i<N; i++) if(D[i]<min_E) min_E = D[i];

	free(Vdag), free(D), free(ham_real);
	return min_E;
}




unsigned long long int choose(int num_electrons)
/*calculating (NUM_SITES choose num_electrons) = NUM_SITES!/((NUM_SITES-num_electrons)!*num_electrons!)*/
{
	int i;
	unsigned long long int c;
	c=1ULL;
	for (i=0;i<num_electrons;i++) c=c*(NUM_SITES-i);
	for (i=0;i<num_electrons;i++) c=c/(i+1);
	return c;
}




int combinations ( int num_electrons,  unsigned long long int *b,int *tab, int N)
/*Returning the combinations of NUM_SITES-choose-elctrons, each represented by NUM_SITES slots of a binary number, where the bit the furthest to the right is site 1, and
furthest to the left is site NUM_SITES. 1 respresents occupied, 0 empty. Finding number of unordered combinations of num_electrons in NUM_SITES.*/
{
	unsigned long long int x,y;
	int i,c,d;
	x=0ULL;
	for (i=0;i<num_electrons;i++) x=x+(1ULL<<i);
	b[0]=x;
	c=0;
	d=0;
	i=0;
	while ((c<NUM_SITES)&& (d<num_electrons))
	{
		if (x & (1ULL<<c))
		{
			tab[i*num_electrons+d]=c+1;
			d++;
		}
		c++;
	}
	for (i=1;i<N;i++)
	{
		y = (x | (x - 1)) + 1;
		x = y | ((((y & -y) / (x & -x)) >> 1) - 1);
		b[i]=x;
		c=0;
		d=0;
		while ((c<NUM_SITES)&& (d<num_electrons))
		{
			if (x & (1ULL<<c))
			{
				tab[i*num_electrons+d]=c+1;
				d++;
			}
			c++;
		}
	}
}




void assign_bonds(int *bonds, int lattice[NX][NY])
/*Attaching a bond number for each bond, allowing assignment of a specific bond-value between two neighbors, where this value is stored in j and k lists*/
{
	int bond_num, site,site2, j, neighbor_count, *neighbors;
	bond_num = 1;
	neighbors = (int*) malloc(4*sizeof(int));
	for(j=0;j<NUM_SITES*NUM_SITES;j++) bonds[j]=0;

	for(site=1; site<NUM_SITES+1;site++)
	{
		neighbor_count = get_neighbors(site, neighbors, lattice, 2);
		for(j=0;j<neighbor_count;j++)
		{
			site2 = neighbors[j];
			if(site2>site)
			{
				bonds[(site-1)*NUM_SITES+site2-1]=bond_num;
				bonds[(site2-1)*NUM_SITES+site-1]=bond_num;
				bond_num++;
			}
		}
	}
	free(neighbors);
}




void copy_arrays(int N, double *k_array, double *j_array, double* b_array,  double* k_to,  double* j_to, double* b_to, int total_steps)
/*storing the updated k, j, and b values*/
{
	memcpy(k_to, k_array, 2*NUM_SITES*total_steps*sizeof(double));
	memcpy(j_to, j_array, 2*NUM_SITES*total_steps*sizeof(double));
	memcpy(b_to, b_array, NUM_SITES*total_steps*sizeof(double));
}



void copy_arrays_bang_bang(double *j_to, double *k_to, double* b_to,  double* j_from,  double* k_from, double* b_from)
{
	memcpy(k_to, k_from, 2*NUMBER_OF_BANGS*sizeof(double));
	memcpy(j_to, j_from, 2*NUMBER_OF_BANGS*sizeof(double));
	memcpy(b_to, b_from, 2*NUMBER_OF_BANGS*sizeof(double));
}


void init_arrays(double *k_array,double *j_array,double *b_array,gsl_rng * r, int total_steps)
/*Initializng the values of the k, j, and b lists which hold the values of the constants for each site (b_array) and between each bond (k_array and j_array)*/
{
	int i,j;
	double random_val1, random_val2, random_val3;
	if(UNIFORM_SITES)
	{
		for (i=0; i<total_steps;i++)
		{

			random_val1 = get_random(LOWJ,UPJ,r);
			random_val2 = get_random(LOWK,UPK,r);
			random_val3 = get_random(LOWB,UPB,r);
			for (j=0; j<NUM_SITES*2; j++)
			{
				j_array[i*NUM_SITES*2+j] = random_val2;
				k_array[i*NUM_SITES*2+j] = random_val3;
			}
			for (j=0; j<NUM_SITES; j++)
			{
				b_array[i*NUM_SITES+j]= random_val1;
			}
		}
	}
	else{
		for (i=0; i<total_steps;i++)
		{
			for (j=0; j<NUM_SITES*2; j++)
			{
				j_array[i*NUM_SITES*2+j] = get_random(LOWJ,UPJ,r);
				k_array[i*NUM_SITES*2+j] = get_random(LOWK,UPK,r);
			}
			for (j=0; j<NUM_SITES; j++)
			{
				b_array[i*NUM_SITES+j]= get_random(LOWB,UPB,r);
			}
		}
	}
	if(PRINTBEST) printf("\nPrinting the best K-J-B arrays from init_arrays"), print_best_arrays(k_array, j_array, b_array, total_steps);
}








void init_arrays_bang(double *j_best_times,double *k_best_times,double *b_best_times,gsl_rng * r)
/*Initializng the values of the k, j, and b lists which hold the values of the constants for each site (b_array) and between each bond (k_array and j_array)*/
{
	int i;
	double random_time_index1, random_time_index2, random_time_index3;
	double tau = TAU_INIT_MCBB;
	k_best_times[0] = 0, j_best_times[0] = 0, b_best_times[0] = 0;
	k_best_times[2*NUMBER_OF_BANGS-1] = tau, j_best_times[2*NUMBER_OF_BANGS-1] = tau, b_best_times[2*NUMBER_OF_BANGS-1] = tau;
	for (i=1; i<NUMBER_OF_BANGS*2-1;i++)
	{

		random_time_index1 = get_random(0,TAU_INIT_MCBB,r);
		random_time_index2 = get_random(0,TAU_INIT_MCBB,r);
		random_time_index3 = get_random(0,TAU_INIT_MCBB,r);
		j_best_times[i] = random_time_index1;
		k_best_times[i] = random_time_index2;
		b_best_times[i] = random_time_index3;
	}
	sort(j_best_times, j_best_times + NUMBER_OF_BANGS*2);
	sort(k_best_times, k_best_times + NUMBER_OF_BANGS*2);
	sort(b_best_times, b_best_times + NUMBER_OF_BANGS*2);

	if(PRINT_TIMES) print_times(j_best_times, k_best_times, b_best_times);
}






void change_array_bang_bang(double *j_times, double *k_times, double *b_times, double change, int random_time_index, int i, double tau)
{
	double upper_bound, lower_bound;
	double *pointer;
	if(i%3 == 0) pointer = j_times;
	if(i%3 == 1) pointer = k_times;
	if(i%3 == 2) pointer = b_times;
//	if(i == 0) printf("TAU FROM CHANGE ARRAY: %f\n", tau);
//	if(i%40 == 0) print_times(j_times, k_times, b_times);
	//print_times(j_times, k_times, b_times);
	//printf("index: %i, Change : %f, i%3: %i\n\n",random_time_index, change, i%3);
	if(change > 0)
	{
		if(random_time_index == 2*NUMBER_OF_BANGS - 1) upper_bound = tau;
		else upper_bound = *(pointer + random_time_index +1);


		if( *(pointer+random_time_index) + change < upper_bound) *(pointer+random_time_index) += change;
		else *(pointer+random_time_index) = upper_bound;
	}
	else
	{
		if(random_time_index == 0) lower_bound = 0;
		else lower_bound = *(pointer+random_time_index-1);


		if( *(pointer+random_time_index) + change > lower_bound) *(pointer+random_time_index) += change;
		else *(pointer+random_time_index) = lower_bound;
	}


/*
	if(i%3 == 1)
	{
		if(change > 0)
		{
			if(j_times[random_time_index] + change <= tau) j_times[random_time_index] += change;
			else j_times[random_time_index] = tau;
		}
		else
		{
			if(j_times[random_time_index] + change >= 0) j_times[random_time_index] += change;
			else j_times[random_time_index] = 0;
		}
		sort(j_times, j_times + NUMBER_OF_BANGS*2);
	}
	if(i%3 == 2)
	{
		if(change > 0)
		{
			if(b_times[random_time_index] + change <= tau) b_times[random_time_index] += change;
			else b_times[random_time_index] = tau;
		}
		else
		{
			if(b_times[random_time_index] + change >= 0) b_times[random_time_index] += change;
			else b_times[random_time_index] = 0;
		}
		sort(b_times, b_times + NUMBER_OF_BANGS*2);
	}
*/
}








void change_array(double *k_array, double *j_array, double *b_array, int random_row, int random_col, double change, int i, int total_steps)
/*changing the j, k, and b arrays. Use the defines at start of program to determine which manipulation functions will be used, and in what order*/
{
	int mod = VARIATIONS;
	if(i%mod==ROW-1) change_row(k_array,j_array,b_array,random_row,change,true, true,true, 1, 0);
	else if(i%mod==COL-1) change_col(total_steps,k_array,j_array,b_array,random_col,change, true, true, true, 1, 0);
	else if(i%mod==ALTROW-1) change_row(k_array,j_array,b_array,random_row,change, true, true, true ,2, 0);
	else if(i%mod==ALTCOL-1) change_col(total_steps,k_array,j_array,b_array,random_col,change, true, true, true, 2, 0);
	else if(i%mod==ALTROW2-1) change_row(k_array,j_array,b_array,random_row,change, true, true, true, 2, 1);
	else if(i%mod==ALTCOL2-1) change_col(total_steps,k_array,j_array,b_array,random_col,change, true, true, true, 2, 1);
	else if(i%mod==ROWK-1) change_row(k_array,j_array,b_array,random_row,change, true, false, false, 1, 0);
	else if(i%mod==ROWJ-1) change_row(k_array,j_array,b_array,random_row,change, false, true, false, 1, 0);
	else if(i%mod==ROWB-1) change_row(k_array,j_array,b_array,random_row,change, false, false, true, 1, 0);
	else if(i%mod==COLK-1) change_col(total_steps,k_array,j_array,b_array,random_col,change, true, false, false, 2, 1);
	else if(i%mod==COLJ-1) change_col(total_steps,k_array,j_array,b_array,random_col,change, false, true, false, 2, 1);
	else if(i%mod==COLB-1) change_col(total_steps,k_array,j_array,b_array,random_col,change, false, false, true, 2, 1);
	else if(i%mod==SINGLE-1) change_single(k_array,j_array,b_array,random_row,random_col,change);
	else printf("NOPE\n");
}



void change_row(double *k_array,double *j_array,double *b_array, int row, double change, bool k, bool j, bool b, int jump, int offset)
/*Changing all of the lists by a value change at at the row number, row. Used in the monte_carlo_simulation_fixed_tau function. offset gives the starting element, jump gives the amount to increase each increment.
  bool k, j, and b determine if the k,j,and b lists will be changes. Bounds of the values that the lists can take are given by the #defines if using a device that's the same as model, +-1 otherwise*/
{
	int i;
	if(k) for (i=offset; i<NUM_SITES*2; i+=jump)
	{
		if ((LOWK < k_array[NUM_SITES*2*row+i] + change) &&  (k_array[NUM_SITES*2*row+i] + change < UPK)) k_array[NUM_SITES*2*row+i] += change;
	}
	if(j) for (i=offset; i<NUM_SITES*2; i+=jump)
	{
		if ((LOWJ < j_array[NUM_SITES*2*row+i] + change) && (j_array[NUM_SITES*2*row+i] + change < UPJ)) j_array[NUM_SITES*2*row+i] += change;
	}
	if(b) for (i=offset; i<NUM_SITES; i+=jump)
	{
		if ((LOWB < b_array[NUM_SITES*row+i] + change) && (b_array[NUM_SITES*row+i] + change < UPB)) b_array[NUM_SITES*row+i] += change;
	}
}




void change_col(int total_steps,double *k_array,double *j_array,double *b_array, int col, double change,bool k, bool j, bool b, int jump, int offset)
/*Changing all of the lists by a value change at at the col number, col. Used in the monte_carlo_simulation_fixed_tau function. offset gives the starting element, jump gives the amount to increase each increment.
  bool k, j, and b determine if the k,j,and b lists will be changes. Bounds of the values that the lists can take are given by the #defines if using a device that's the same as model, +-1 otherwise*/
{
	int i;

	if(k) for (i=offset; i<total_steps; i+=jump)
	{
		if ((LOWK < k_array[NUM_SITES*2*i+col] + change) && (k_array[NUM_SITES*2*i+col] + change  < UPK)) k_array[NUM_SITES*2*i+col] += change;
	}
	if(j) for (i=offset; i<total_steps; i+=jump)
	{
		if ((LOWJ < j_array[NUM_SITES*2*i+col] + change) && (j_array[NUM_SITES*2*i+col] + change < UPJ)) j_array[NUM_SITES*2*i+col] += change;
	}
	if(b) for (i=offset; i<total_steps; i+=jump)
	{
		if ((LOWB < b_array[NUM_SITES*i+(int)floor(col/2.0)] + change) && (b_array[NUM_SITES*i+(int)floor(col/2.0)] + change < UPB)) b_array[NUM_SITES*i+(int)floor(col/2.0)] += change;
	}
}




void change_single(double *k_array,double *j_array,double *b_array, int row,int col, double change)
/*Changing a single element in each array, at column col, row row*/
{
	if ((LOWJ < j_array[NUM_SITES*2*row+col] + change) && (j_array[NUM_SITES*2*row+col] + change < UPJ)) j_array[NUM_SITES*2*row+col] += change;
 	if ((LOWK < k_array[NUM_SITES*2*row+col] + change) && (k_array[NUM_SITES*2*row+col] + change < UPK)) k_array[NUM_SITES*2*row+col] += change;
	if ((LOWB < b_array[NUM_SITES*row+(int)floor(col/2.0)] + change) && (b_array[NUM_SITES*row+(int)floor(col/2.0)] + change  < UPB)) b_array[NUM_SITES*row+(int)floor(col/2.0)] += change;
}



void scale_arrays(int total_steps, double* j_best, double* k_best,double* b_best)
{
	int i, j;
	for(i=total_steps;i>0;i--)
	{
		for(j=0;j<NUM_SITES*2;j++)
		{

			k_best[NUM_SITES*2*(ARRAY_SCALAR*i-1)+j] = 	k_best[NUM_SITES*2*(i-1)+j];
			k_best[NUM_SITES*2*(ARRAY_SCALAR*i-2)+j] = 	k_best[NUM_SITES*2*(i-1)+j];
			j_best[NUM_SITES*2*(ARRAY_SCALAR*i-1)+j] = 	j_best[NUM_SITES*2*(i-1)+j];
			j_best[NUM_SITES*2*(ARRAY_SCALAR*i-2)+j] = 	j_best[NUM_SITES*2*(i-1)+j];

		}
		for(j=0;j<NUM_SITES;j++)
		{

			b_best[NUM_SITES*(ARRAY_SCALAR*i-1)+j] = 	b_best[NUM_SITES*(i-1)+j];
			b_best[NUM_SITES*(ARRAY_SCALAR*i-2)+j] = 	b_best[NUM_SITES*(i-1)+j];

		}
	}
}





void construct_lattice(int lattice[NX][NY])
/*Build a NX*NX lattice, with sites 1 through NX*NX listed in a snaking pattern*/
{
	int x,y;
	for (x=0;x<NX;x++)
	{
		if (x%2 ==1) for (y=0;y<NY;y++) lattice[x][y] = (NX*x)+y+1;
		else for (y=0;y<NY;y++) lattice[x][NY-1-y] = NX*x+y+1;
	}
}




int get_neighbors(int site, int *neighbors, int lattice[NX][NY], int D)
/*Gets the neighbors of a site, returning all 4 neighbors if open_boundry coditions are true, otherwise just closed-boudnary neighbors*/
{
	int x,y,i,count=0;
	if (D==2)
	{
		for (x=0;x<NX;x++) for (y=0;y<NY;y++) if(site == lattice[x][y]) goto end_loop;//Finding the coordinates
		end_loop: for(i=0;i<4;i++) neighbors[i]=0;

		if (PERIODIC)
		{
			neighbors[0] = lattice[(x+1)%NX][y];
			neighbors[1] = lattice[(x+NX-1)%NX][y];
			neighbors[2] = lattice[x][(y+1)%NY];
			neighbors[3] = lattice[x][(y+(NY-1))%NY];
			count = 4;
		}
		else
		{
			if (x+1<NX) neighbors[count] = lattice[x+1][y], count++;
			if (x>0) neighbors[count] = lattice[x-1][y], count++;
			if (y+1<NY) neighbors[count] = lattice[x][y+1], count++;
			if (y>0) neighbors[count] = lattice[x][y-1], count++;
		}
		return count;
	}
	else if (D==1)
	{
		if (PERIODIC)
		{
			neighbors[0] = (site%(NX*NY))+1;
			neighbors[1] = ((site+(NX*NY)-2)%(NX*NY))+1;
			count = 2;
		}
		else
		{
			if (site<(NX*NY)) neighbors[count] = site+1, count++;
			if (site>1) neighbors[count] = site-1, count++;
		}
		return count;
	}
	printf("\n\n\nERROR! NO NEIGHBORS FOUND.\n\n\n");
}




double get_random(double lower, double upper, gsl_rng * r)
/*randomly generating a number in [lower, upper), including lower, up to upper (but not including)*/
{
	double u;

	/*varying seed, based on time
	const gsl_rng_type * TT;
	gsl_rng * r;
	struct timeval tv;
	double seed;
	gsl_rng_env_setup();
	gettimeofday(&tv,0);
	seed = tv.tv_usec;
	TT = gsl_rng_default;
	r = gsl_rng_alloc (TT);
	gsl_rng_set(r, seed);
	u = gsl_rng_uniform(r);
	gsl_rng_free (r);
	return u*(upper-lower)+lower;
	*/

	u = gsl_rng_uniform(r);
	return u*(upper-lower)+lower;
}





void check_norm(double* psi, int N)
{
	int i,j=0;
	double sum=0;
	for(i=0;i<N*2; i+=2)
	{
		sum+= psi[i]*psi[i]+(psi[i+1]*psi[i+1]);

	}
	if(sum>1.000000001 or sum<0.99999999) printf("\n\n\nNORM ERROR, SIZE: %f\n\n\n", sum);
}




void check_unitary(double* ham, int N)
{
	int i,j, *IPIV, LWORK=N*N, INFO, LDA=N, LDB=N, LDC=N;
	double *ham_t,*unitary, *WORK,ALPHA[2], BETA[2];
	char TRANSA = 'C', TRANSB = 'N';


	ALPHA[0]=1.0, ALPHA[1]=0.0;
	BETA[0]=0.0, BETA[1]=0.0;

	ham_t = (double*) malloc(2*N*N*sizeof(double));
	unitary = (double*) malloc(2*N*N*sizeof(double));
	memcpy(ham_t, ham, sizeof(double)*2*N*N);
	zgemm_(&TRANSA, &TRANSB, &N, &N, &N, ALPHA, ham_t, &LDA, ham, &LDB, BETA, unitary, &LDC); //matrix mult

	for(i=0;i<N;i++)
	{
		unitary[2*(i*N+i)] = unitary[2*(i*N+i)] -1;
	}

	for(i=0;i<N*N*2;i++)
	{
		if(unitary[i] < -0.00000001 or unitary[i] > 0.00000001) printf("\n\n\nERROR, NON UNITARY ELEMENTS AT %i, VALUE: %f\n\n\n", i, unitary[i]);

	}
	free(ham_t), free(unitary);
}





void check_hermicity(double* ham, int N)
{
	int i,j;
	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{
			if(abs(ham[2*(j*N+i)] - ham[2*(i*N+j)]) > 0.00001) printf("\n\n\nERROR, NON HERMITIAN ELEMENT AT i,j: %i,%i\nElement_ij = %10.7f\nElement_ji = %10.7f\n\n\n", i*N+j, j*N+i, ham[2*(i*N+j)], ham[2*(j*N+i)]);
		}
	}
}





void check_weights(double* state, double* ham, int N)
{
	int i,j;
	double *D, *Vdag, *ham_real;
	Vdag = (double*) malloc(N*N*sizeof(double));
	D = (double*) malloc(N*sizeof(double));
	ham_real = (double*) malloc(N*N*sizeof(double));
	for (i=0;i<N*N;i++) ham_real[i] = ham[2*i];

	diag_hermitian_real_double(N, ham_real,Vdag, D);

	double sum_real, sum_im, c_squared_sum=0;

	for(j=0;j<N;j++)
	{
		sum_real=0;
		sum_im=0;
		for(i=0; i<N; i++)
		{
			sum_real += state[2*i]*Vdag[j*N+i];
			sum_im += -state[2*i+1]*Vdag[j*N+i];
		}
		c_squared_sum += sum_real*sum_real+sum_im*sum_im;
	}
	if(c_squared_sum>1.00001 or c_squared_sum <0.999999) printf("\n\n\nERROR, BAD WEIGHTS\n\n\n");
	free(Vdag), free(D), free(ham_real);
}



void export_mc_data(double pre_converge_tau, int pre_converge_steps, double pre_converge_E, double *j_best, double *k_best, double *b_best,double *pre_converge_j_best, double *pre_converge_k_best, double *pre_converge_b_best, double *tau_array, double *best_E_array, int total_steps, double tau, double ground_E, double *jkb_initial, double *jkb_target, int index, int seed, int num_electrons, double *gf)
{
	int i;
	ofstream file;

	string uni, pbc, gi, gt, fi, ft;

	gi = to_string(gf[0]);
        gi.erase(gi.find_last_not_of('0') + 2, string::npos);
	fi = to_string(gf[1]);
        fi.erase(fi.find_last_not_of('0') + 2, string::npos);
	gt = to_string(gf[2]);
        gt.erase(gt.find_last_not_of('0') + 2, string::npos);
	ft = to_string(gf[3]);
        ft.erase(ft.find_last_not_of('0') + 2, string::npos);

	if (PERIODIC) pbc = 't';
	else pbc = 'f';
	if (UNIFORM_SITES) uni = 't';
	else uni = 'f';

	string dir = "data/" + to_string(NX) + "x" +to_string(NY) + "/" + to_string(num_electrons) ;
	string file_name = "_occupants/MC_____PBC="+pbc+"_UNI="+uni+"_DD="+to_string(DEVICE_DIMENSION)+"___gi=" + gi + "_gt=" + gt + "_fi="+ fi + "_ft=" + ft +".txt";

	if(seed == 1)
	{
		file.open(dir + file_name);
		file << "##############################################################\n";
		file << "SIMULATIONS PARAMETERS\n";
		file << "##############################################################\n";
		file << "DIAG=                       " <<  boolalpha << DIAG << "\n";
		file << "SEED_TOTAL=                 " <<  SEED_TOTAL << "\n";
		file << "DIFFERENCE_LIMIT_MC=        " <<  DIFFERENCE_LIMIT_MC << "\n";
		file << "TAU_INIT_MC=                " <<  TAU_INIT_MC << "\n";
		file << "MAX_TAU_MC=                 " <<  MAX_TAU_MC << "\n";
		file << "TAU_SCALAR_MC=              " <<  TAU_SCALAR_MC << "\n";
		file << "MAX_CHANGE_MC_INIT=         " <<  MAX_CHANGE_MC_INIT << "\n";
		file << "ACCEPTANCE_PROB_MC=         " <<  ACCEPTANCE_PROB_MC << "\n";
		file << "TEMP_EXP_DECAY_MC=          " <<  TEMP_EXP_DECAY_MC << "\n";
		file << "BINARY_SEARCH_TAU_LIMIT_MC= " <<  BINARY_SEARCH_TAU_LIMIT_MC << "\n";
		file << "RANDOM_STATES_MC=           " <<  RANDOM_STATES_MC  << "\n";
		file << "SWEEPS_MC=                  " <<  SWEEPS_MC << "\n";
		file << "TOTAL_STEPS_INIT_MC=        " <<  TOTAL_STEPS_INIT_MC << "\n";
		file << "TEMP_DECAY_ITERATIONS_MC=   " <<  TEMP_DECAY_ITERATIONS_MC  << "\n";
		file << "TEMP_DECAY_LIMIT_MC=        " <<  TEMP_DECAY_LIMIT_MC << "\n";
		file << "MAX_EVOLVE_STEPS_MC=        " <<  MAX_EVOLVE_STEPS_MC << "\n";
		file << "MAX_TAU_STEPS_MC=           " <<  MAX_TAU_STEPS_MC << "\n";
		file << "ARRAY_SCALAR=               " <<  ARRAY_SCALAR << "\n\n\n";

		file << "##############################################################\n";
		file << "HAMILTONIAN PARAMETERS\n";
		file << "##############################################################\n";
		file << "g_initial=        " << gi << "\n";
		file << "f_initial=        " << fi << "\n";
		file << "g_target=         " << gt << "\n";
		file << "f_target=         " << ft << "\n";
		file << "j_initial=        " << jkb_initial[0] << "\n";
		file << "k_initial=        " << jkb_initial[1] << "\n";
		file << "b_initial=        " << jkb_initial[2] << "\n";
		file << "j_target=         " << jkb_target[0] << "\n";
		file << "k_target=         " << jkb_target[1] << "\n";
		file << "b_target=         " << jkb_target[2] << "\n";
		file << "PERIODIC=         " << boolalpha << PERIODIC << "\n";
		file << "UNIFORM_SITES=    " << boolalpha << UNIFORM_SITES << "\n";
		file << "DEVICE_DIMENSION= " << DEVICE_DIMENSION << "\n";
		file << "MAX_PARAM=        " << MAX_PARAM << "\n";
		file << "MIN_PARAM=        " << MIN_PARAM << "\n\n";
	}
	else file.open(dir + file_name, ios::app);


	file << "\n\n\n\n##############################################################################\n";
	file << "##############################################################################\n";
	file << "seed =                     " << seed << "\n";
	file << "tau =                      " << tau   << "\n";
	file << "total_steps =              " << total_steps << "\n\n";

	file << "pre_converge_tau =         " << pre_converge_tau << "\n";
	file << "pre_converge_E =           " << pre_converge_E << "\n";
	file << "pre_converge_total_steps = " << pre_converge_steps << "\n";
	file << "pre_converge_j_protocol = [";
	for(i=0;i<pre_converge_steps;i++) file << pre_converge_j_best[i*NUM_SITES*2] << ", ";
	file << "]\npre_converge_k_protocol = [";
	for(i=0;i<pre_converge_steps;i++) file << pre_converge_k_best[i*NUM_SITES*2] << ", ";
	file << "]\npre_converge_b_protocol = [";
	for(i=0;i<pre_converge_steps;i++) file << pre_converge_b_best[i*NUM_SITES] << ", ";
	file << "]\n\n";

	file << "j_protocol = [";
	for(i=0;i<total_steps;i++) file << j_best[i*NUM_SITES*2] << ", ";
	file <<"]\nk_protocol = [";
	for(i=0;i<total_steps;i++) file << k_best[i*NUM_SITES*2] << ", ";
	file <<"]\nb_protocol = [";
	for(i=0;i<total_steps;i++) file << b_best[i*NUM_SITES] << ", ";
	file <<"]\n\n";

	file <<    "times =    [";
	for(i=0;i<index+1;i++) file << tau_array[i] << ", ";
	file << "]\nE_best =   [";
	for(i=0;i<index+1;i++) file << best_E_array[i] << ", ";
	file << "]\nE_ground = [";
	for(i=0;i<index+1;i++) file << ground_E << ", ";
	file << "]\n\n";

	file.close();
}








void export_mcbb_data(double* pre_converge_taus,double* pre_converge_Es,double* pre_converge_j_times,double*  pre_converge_k_times,double* pre_converge_b_times, double *j_best_times, double *k_best_times, double *b_best_times, double *tau_array, double *best_E_array, double tau, double ground_E, double *jkb_initial, double *jkb_target, int index, int seed, int num_electrons, double *gf)
{
	int i;
	ofstream file;

	string uni, pbc, gi, gt, fi, ft;

	gi = to_string(gf[0]);
        gi.erase(gi.find_last_not_of('0') + 2, string::npos);
	fi = to_string(gf[1]);
        fi.erase(fi.find_last_not_of('0') + 2, string::npos);
	gt = to_string(gf[2]);
        gt.erase(gt.find_last_not_of('0') + 2, string::npos);
	ft = to_string(gf[3]);
        ft.erase(ft.find_last_not_of('0') + 2, string::npos);




	if (PERIODIC) pbc = 't';
	else pbc = 'f';
	if (UNIFORM_SITES) uni = 't';
	else uni = 'f';

	string dir = "data/" + to_string(NX) + "x" +to_string(NY) + "/" + to_string(num_electrons) ;
	string file_name = "_occupants/MCBB___PBC="+pbc+"_UNI="+uni+"_DD="+to_string(DEVICE_DIMENSION)+"___gi=" + gi + "_gt=" + gt + "_fi="+ fi + "_ft="+ ft +".txt";

	if(seed == 1)
	{
		file.open(dir + file_name);
		file << "##############################################################\n";
		file << "SIMULATIONS PARAMETERS\n";
		file << "##############################################################\n";
		file << "DIAG=                         " <<  boolalpha << DIAG << "\n";
		file << "SEED_TOTAL=                   " <<  SEED_TOTAL << "\n";
		file << "DIFFERENCE_LIMIT_MCBB=        " <<  DIFFERENCE_LIMIT_MCBB << "\n";
		file << "TAU_INIT_MCBB=                " <<  TAU_INIT_MCBB << "\n";
		file << "MAX_TAU_MCBB=                 " <<  MAX_TAU_MCBB << "\n";
		file << "TAU_SCALAR_MCBB=              " <<  TAU_SCALAR_MCBB << "\n";
		file << "CHANGE_FRACTION_MCBB=         " <<  CHANGE_FRACTION_MCBB << "\n";
		file << "ACCEPTANCE_PROB_MCBB=         " <<  ACCEPTANCE_PROB_MCBB << "\n";
		file << "TEMP_EXP_DECAY_MCBB=          " <<  TEMP_EXP_DECAY_MCBB << "\n";
		file << "BINARY_SEARCH_TAU_LIMIT_MCBB= " <<  BINARY_SEARCH_TAU_LIMIT_MCBB << "\n";
		file << "RANDOM_STATES_MCBB=           " <<  RANDOM_STATES_MCBB  << "\n";
		file << "NUMBER_OF_BANGS=              " <<  NUMBER_OF_BANGS << "\n";
		file << "SWEEPS_MCBB=                  " <<  SWEEPS_MCBB << "\n";
		file << "TEMP_DECAY_ITERATIONS_MCBB=   " <<  TEMP_DECAY_ITERATIONS_MCBB  << "\n";
		file << "TEMP_DECAY_LIMIT_MCBB=        " <<  TEMP_DECAY_LIMIT_MCBB << "\n";
		file << "MAX_TAU_STEPS_MCBB=           " <<  MAX_TAU_STEPS_MCBB << "\n";
		file << "CONVERGE_LIMITS=             [";
		for(i=0;i<sizeof(CONVERGE_LIMITS)/sizeof(double);i++) file << CONVERGE_LIMITS[i] << ", ";
		file << "]\n\n";




		file << "##############################################################\n";
		file << "HAMILTONIAN PARAMETERS\n";
		file << "##############################################################\n";
		file << "g_initial=        " << gi << "\n";
		file << "f_initial=        " << fi << "\n";
		file << "g_target=         " << gt << "\n";
		file << "f_target=         " << ft << "\n";
		file << "j_initial=        " << jkb_initial[0] << "\n";
		file << "k_initial=        " << jkb_initial[1] << "\n";
		file << "b_initial=        " << jkb_initial[2] << "\n";
		file << "j_target=         " << jkb_target[0] << "\n";
		file << "k_target=         " << jkb_target[1] << "\n";
		file << "b_target=         " << jkb_target[2] << "\n";
		file << "PERIODIC=         " << boolalpha << PERIODIC << "\n";
		file << "UNIFORM_SITES=    " << boolalpha << UNIFORM_SITES << "\n";
		file << "DEVICE_DIMENSION= " << DEVICE_DIMENSION << "\n";
		file << "MAX_PARAM=        " << MAX_PARAM << "\n";
		file << "MIN_PARAM=        " << MIN_PARAM << "\n\n";
	}
	else file.open(dir + file_name, ios::app);


	file << "\n\n\n\n##############################################################################\n";
	file << "##############################################################################\n";
	file << "seed =             " << seed << "\n";
	file << "tau =              " << tau   << "\n";

	file << "\ninitial_tau=      " << pre_converge_taus[0];
	file << "\ninitial_E=        " << pre_converge_Es[0];
	file << "\ninitial_j_times = [";
	for(i=0;i<2*NUMBER_OF_BANGS;i++) file << pre_converge_j_times[i] << ", ";
	file << "]\ninitial_k_times = [";
	for(i=0;i<2*NUMBER_OF_BANGS;i++) file << pre_converge_k_times[i] << ", ";
	file << "]\ninitial_b_times = [";
	for(i=0;i<2*NUMBER_OF_BANGS;i++) file << pre_converge_b_times[i] << ", ";
	file << "]\n";

	file << "\npre_converge_taus = [";
	for(i=1;i<sizeof(CONVERGE_LIMITS)/sizeof(double);i++) file << pre_converge_taus[i] << ", ";
	file << "]\npre_converge_Es =  [";
	for(i=1;i<sizeof(CONVERGE_LIMITS)/sizeof(double);i++) file << pre_converge_Es[i] << ", ";
	file << "]\npre_converge_j_times = [";
	for(i=2*NUMBER_OF_BANGS;i<2*NUMBER_OF_BANGS*sizeof(CONVERGE_LIMITS)/sizeof(double);i++) file << pre_converge_j_times[i] << ", ";
	file << "]\npre_converge_k_times = [";
	for(i=2*NUMBER_OF_BANGS;i<2*NUMBER_OF_BANGS*sizeof(CONVERGE_LIMITS)/sizeof(double);i++) file << pre_converge_k_times[i] << ", ";
	file << "]\npre_converge_b_times = [";
	for(i=2*NUMBER_OF_BANGS;i<2*NUMBER_OF_BANGS*sizeof(CONVERGE_LIMITS)/sizeof(double);i++) file << pre_converge_b_times[i] << ", ";
	file << "]\n";




	file << "\nj_times = [";
	for(i=0;i<2*NUMBER_OF_BANGS;i++) file << j_best_times[i] << ", ";
	file << "]\nk_times = [";
	for(i=0;i<2*NUMBER_OF_BANGS;i++) file << k_best_times[i] << ", ";
	file << "]\nb_times = [";
	for(i=0;i<2*NUMBER_OF_BANGS;i++) file << b_best_times[i] << ", ";
	file << "]\nj/k/b = [";
	for(i=1;i<2*NUMBER_OF_BANGS+1;i++) file << i%2 << ", ";
	file <<"]\n\n";

	file <<    "times =    [";
	for(i=0;i<index+1;i++) file << tau_array[i] << ", ";
	file << "]\nE_best =   [";
	for(i=0;i<index+1;i++) file << best_E_array[i] << ", ";
	file << "]\nE_ground = [";
	for(i=0;i<index+1;i++) file << ground_E << ", ";
	file << "]\n";

	file.close();
}





void export_adiabatic_data(double *tau_array, double *E_array, int total_steps, double tau, double ground_E, double g_initial, double f_initial, double g_target, double f_target, double *jkb_initial, double*jkb_target, int index, int num_electrons)
{
	int i;
	ofstream file;


	string uni, pbc, gi, gt, f;


	gi = to_string(g_initial);
        gi.erase(gi.find_last_not_of('0') + 2, string::npos);
	gt = to_string(g_target);
        gt.erase(gt.find_last_not_of('0') + 2, string::npos);
	f = to_string(f_initial);
        f.erase(f.find_last_not_of('0') + 2, string::npos);


	if (PERIODIC) pbc = 't';
	else pbc = 'f';
	if (UNIFORM_SITES) uni = 't';
	else uni = 'f';


	string dir = "data/" + to_string(NX) + "x" +to_string(NY) + "/" + to_string(num_electrons) ;
	string file_name = "_occupants/ADIA___PBC="+pbc+"_UNI="+uni+"_DD="+to_string(DEVICE_DIMENSION)+"___gi=" + gi + "_gt=" + gt + "_f="+ f +".txt";


	file.open(dir + file_name);
	file << "##############################################################\n";
	file << "SIMULATIONS PARAMETERS\n";
	file << "##############################################################\n";
	file << "DIAG=                       " <<  boolalpha << DIAG << "\n";
	file << "DIFFERENCE_LIMIT_ADIA=      " <<  DIFFERENCE_LIMIT_ADIA << "\n";
	file << "TAU_INIT_ADIA=              " <<  TAU_INIT_ADIA << "\n";
	file << "MAX_TAU_ADIA=               " <<  MAX_TAU_ADIA << "\n";
	file << "TAU_SCALAR_ADIA=            " <<  TAU_SCALAR_ADIA << "\n";
	file << "TIME_STEP_ADIA=             " <<  TIME_STEP_ADIA << "\n\n";


	file << "##############################################################\n";
	file << "HAMILTONIAN PARAMETERS\n";
	file << "##############################################################\n";
	file << "g_initial=        " << gi << "\n";
	file << "f_initial=        " << f << "\n";
	file << "g_target=         " << gt << "\n";
	file << "f_target=         " << f << "\n";
	file << "j_initial=        " << jkb_initial[0] << "\n";
	file << "k_initial=        " << jkb_initial[1] << "\n";
	file << "b_initial=        " << jkb_initial[2] << "\n";
	file << "j_target=         " << jkb_target[0] << "\n";
	file << "k_target=         " << jkb_target[1] << "\n";
	file << "b_target=         " << jkb_target[2] << "\n";
	file << "PERIODIC=         " << boolalpha << PERIODIC << "\n";
	file << "UNIFORM_SITES=    " << boolalpha << UNIFORM_SITES << "\n";
	file << "DEVICE_DIMENSION= " << DEVICE_DIMENSION << "\n";
	file << "MAX_PARAM=        " << MAX_PARAM << "\n";
	file << "MIN_PARAM=        " << MIN_PARAM << "\n\n";




	file << "\n\n\n\n##############################################################################\n";
	file << "##############################################################################\n";
	file << "tau =              " << tau   << "\n";
	file << "total_steps " << total_steps << "\n\n";

	file << "times= [";
	for(i=0;i<index;i++) file << tau_array[i] << ", ";
	file << "]\nE_best = [";
	for(i=0;i<index;i++) file << E_array[i] << ", ";
	file << "]\nE_ground = [";
	for(i=0;i<index;i++) file << ground_E << ", ";
	file << "]";
	file.close();
}






void export_evolve_data(int *table, unsigned long long int *b, int num_electrons, int N,double* j_best, double* k_best, double* b_best, double tau, double time_step, int total_steps, int lattice[NX][NY],double* psi_start, double ground_E, double* ham_target, int seed)
{
	ofstream file;
	file.open("data/exp_difference_vs_time/E_VS_T___ToteTime=" + to_string(tau) + "_ToteStep" + to_string(total_steps) + "_SEED=" + to_string(seed) + ".txt");
	char output[200];
	sprintf(output, "TOTAL TIME  |  TIME STEP  |  TOTAL STEPS\n %5f      |%10f   |   %5i\n", tau, time_step,total_steps);
	file << output;

	int *bonds, i, j;
	double *k_avg, *b_avg, *j_avg, *E_list, *T_list, *psi, j_av,k_av,b_av;

	bonds = (int*) malloc(NUM_SITES*NUM_SITES*sizeof(int));
	E_list = (double*) malloc(total_steps*sizeof(double));
	T_list = (double*) malloc(total_steps*sizeof(double));
	psi = (double*) malloc(N*2*sizeof(double));
	k_avg = (double*) malloc(total_steps*sizeof(double));
	j_avg = (double*) malloc(total_steps*sizeof(double));
	b_avg = (double*) malloc(total_steps*sizeof(double));

	assign_bonds(bonds, lattice);
	memcpy(psi,psi_start, 2*N*sizeof(double));


	for(i=0;i<total_steps;i++) //Averaging all of the j-K-B values
	{
		k_av =0, b_av =0, j_av=0;

		for(j=0;j<NUM_SITES*2;j++) j_av += j_best[i*NUM_SITES*2+j], k_av += k_best[i*NUM_SITES*2+j];
		for(j=0;j<NUM_SITES;j++) b_av += b_best[i*NUM_SITES+j];

		k_avg[i] = k_av/(NUM_SITES*2);
		j_avg[i] = j_av/(NUM_SITES*2);
		b_avg[i] = b_av/(NUM_SITES);
	}

	file << "\n\n\nK averages\n[";
	for(i=0;i<total_steps;i++) file << k_avg[i] << ", ";

	file << "]\n\n\nJ averages\n[";
	for(i=0;i<total_steps;i++) file << j_avg[i] << ", ";

	file << "]\n\n\nB averages\n[";
	for(i=0;i<total_steps;i++) file << b_avg[i] << ", ";
	file << "]\n\n\n";

	export_evolve_calc(table,b,num_electrons,N,lattice, psi, k_best,j_best,b_best, bonds, total_steps, time_step, E_list, T_list, ham_target);

	file << "\n\n\ntimes=\n[";
	for(i=0;i<total_steps;i++) file << T_list[i] << ", ";

	file << "]\n\n\nExpectation-E_ground=\n[";
	for(i=0;i<total_steps;i++) file << (E_list[i]-ground_E) << ", ";
	file << "]";


	free(k_avg), free(b_avg), free(j_avg), free(T_list), free(E_list), free(bonds), free(psi);
	file.close();
}



void export_evolve_calc(int *table, unsigned long long int *b,int num_electrons,int N,  int lattice[NX][NY], double *psi, double *k_array, double *j_array, double *b_array, int *bonds, int total_steps, double time_step, double* E_list, double* T_list, double* ham_target)
{
	int i,j;
	double *ham_dev,*ham_t_i, *ham_real,*exp_matrix,*D, *Vdag;
	ham_dev = (double *) malloc (2*N*N*sizeof (double));
	exp_matrix = (double *) malloc (2*N*N*sizeof (double));
	ham_t_i = (double *) malloc (2*N*N*sizeof (double));
	ham_real = (double *) malloc(N*N*sizeof(double));
	Vdag = (double*) malloc(N*N*sizeof(double));
	D = (double*) malloc(N*sizeof(double));

	for (i=0; i<total_steps;i++)
	{
		construct_device_hamiltonian(table, b, num_electrons, N, ham_dev, lattice, k_array, j_array, b_array, bonds,i,DEVICE_DIMENSION);

		if(CHECK) check_norm(psi, N);

		if(DIAG)
		{
			for (j=0; j<N*N; j++) ham_real[j]=0.0,Vdag[j]=0.0;
			for (j=0; j<N; j++) D[j]=0.0;
			for (j=0; j<N*N; j++) ham_real[j] = ham_dev[2*j];//converting an all-real-valued complex matrix into just real matrix
			diag_hermitian_real_double(N, ham_real,Vdag, D);
			exp_diaganolized_mat(ham_dev, Vdag, D, N, time_step);//This function exponentiates D to e^(-i*time_step*D)

			if(CHECK) check_unitary(ham_dev, N);
			matrix_vector_mult(ham_dev,psi, N);
		}
		else
		{
			for (j=0; j<N*N*2; j++) exp_matrix[j] = 0.0, ham_t_i[j]=0.0;
			for (j=0; j<N*N; j++)
			{
				ham_t_i[2*j+1] = (ham_dev[2*j]*-time_step); //multiplying by -i*dt for the Pade approximation
				ham_t_i[2*j] = (ham_dev[2*j+1]*time_step);
			}
			exp_general_complex_double(N, ham_t_i, exp_matrix);
			if(CHECK) check_unitary(exp_matrix, N);
			matrix_vector_mult(exp_matrix, psi, N);
		}
		if(CHECK) check_norm(psi, N);
		E_list[i] = cost(psi, ham_target, N);
		T_list[i] = (i+1)*time_step;

	}
	free(ham_dev), free(ham_t_i), free(exp_matrix), free(D), free(Vdag),free(ham_real);
}




void print_vector(double* psi,int N)
{
	int i;
	printf("[");
	for(i=0;i<N;i++)
	{
		printf("%04.3f+%04.3fi; ", psi[2*i], psi[2*i+1]);
	}
	printf("]\n");
}




void print_best_arrays(double *k_array, double *j_array, double *b_array, int total_steps)
{
	int i,j;

	if (UNIFORM_SITES)
	{
		printf("\nK_array= %iX%i (stepsXsites)\n[", total_steps, NUM_SITES*2);
		for(i=0;i<total_steps;i++) printf("%5.2f ",k_array[i*NUM_SITES*2]);
		printf(";]\n");

		printf("J_array= %iX%i (stepsXsites)\n[", total_steps, NUM_SITES*2);
		for(i=0;i<total_steps;i++) printf("%5.2f ",j_array[i*NUM_SITES*2]);
		printf(";]\n");

		printf("B_array= %iX%i (stepsXsites)\n[", total_steps, NUM_SITES);
		for(i=0;i<total_steps;i++) printf("%5.2f ",b_array[i*NUM_SITES]);
		printf(";]\n");
	}

	else
	{
		printf("\nK_array= %iX%i (stepsXsites)\n[", total_steps, NUM_SITES*2);
		for(i=0;i<total_steps;i++)
		{
			for(j=0;j<2*NUM_SITES;j++)
			{
				printf("%5.2f ",k_array[i*NUM_SITES*2+j]);
				if(k_array[i*NUM_SITES+j] > 1 || -1 > k_array[i*NUM_SITES+j]) printf("\n\n\nERROR IN PRINT BEST\n\n\n");
			}
			if (i==total_steps-1) printf("]\n");
			else printf(";\n ");
		}

		printf("J_array= %iX%i (stepsXsites)\n[", total_steps, NUM_SITES*2);
		for(i=0;i<total_steps;i++)
		{
			for(j=0;j<2*NUM_SITES;j++)
			{
				printf("%5.2f ",j_array[i*NUM_SITES*2+j]);
				if(k_array[i*NUM_SITES+j] > 1 || -1 > k_array[i*NUM_SITES+j]) printf("\n\n\nERROR IN PRINT BEST\n\n\n");
			}
			if (i==total_steps-1) printf("]\n");
			else printf(";\n ");
		}

		printf("B_array= %iX%i (stepsXsites)\n[", total_steps, NUM_SITES);
		for(i=0;i<total_steps;i++)
		{
			for(j=0;j<NUM_SITES;j++)
			{
				printf("%5.2f ",b_array[i*NUM_SITES+j]);
				if(b_array[i*NUM_SITES+j] > 1 || -1 > b_array[i*NUM_SITES+j]) printf("\n\n\nERROR IN PRINT BEST\n\n\n");
			}
			if (i==total_steps-1) printf("]\n");
			else printf(";\n ");
		}
		printf("\n");
	}
}





void print_hamiltonianR(double* hamiltonian, int N)
/*prints a real-valued hamiltonian in matrix form*/
{
	int i,j;
	printf("\nPrinting the %ix%i Hamiltonian=\n[", N,N);
	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++) printf("%09.5f ",(hamiltonian[j*N+i])+0.0);
		if (i==N-1) printf("]");
		else printf("\n ");
	}
	printf("\n");
}




void print_hamiltonian(double* hamiltonian, int N)
/*prints a complex hamiltonian in matrix form, with the option to just print out the real values*/
{
	int i,j;
	printf("\nPrinting the %ix%i Hamiltonian=\n[", N,N);
	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++) printf("%09.5f+%09.5fi  ",(hamiltonian[2*(j*N+i)]+0.0), hamiltonian[2*(j*N+i)+1]);
		if (i==N-1) printf("]");
		else printf(";\n ");
	}
	printf("\n");
}




void print_E(double ground_E, double best_E, double initial_E)
{
	printf("Post-Monte_Carlo Expectation:  %f\n", best_E);
	printf("#######################################################################\n");
	printf("BEST ENERGY=    %9.6f\n",best_E);
	printf("TARGET ENERGY=  %9.6f\n",ground_E);
	printf("DIFFERENCE=     %9.6f\n",best_E - ground_E);
	printf("INITIAL ENERGY= %9.6f\n",initial_E);
	printf("DISTANCNCE=     %9.6f\n",(best_E - ground_E) / (initial_E - ground_E));
	printf("#######################################################################\n");
}




void test_function()
{
	int i,j;
	int N = 3;
	double *D, *Vdag, *ham_target, min_E=1000;
	Vdag = (double*) malloc(N*N*sizeof(double));
	D = (double*) malloc(N*sizeof(double));
}
