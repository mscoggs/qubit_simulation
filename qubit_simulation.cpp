#include <iomanip>
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

//g++ -o qubit_simulation qubit_simulation.cpp -llapack -lblas -lgsl
//./qubit_simulation

//Physical Parameters
#define OPEN true
#define DEVICE_DIMENSION 2
#define T 1
#define V 2
#define J_START 0
#define K_START 0
#define B_START 1


//Non-Physical Parameters
#define SEED 1
#define NX 2
#define NY NX
#define NUM_SITES NX*NY
#define DEVICE_SAME_MODEL true
#define DIAG true
#define LIMIT 30
#define RANDOM_STATES 3
#define SWEEPS 200
#define CHANGE 0.02
#define ACCEPTANCE_PROB 0.7
#define EXP_DECAY 0.90
#define TEMP_DECAY_ITERATIONS 5
#define UNIFORM true

#ifdef DEVICE_SAME_MODEL
#define UPJ 1//max(J_START, J_TARGET)
#define UPK 1//max(K_START, K_TARGET)
#define UPB 1//max(B_START, B_TARGET)
#define LOWJ 0//min(J_START, J_TARGET)
#define LOWK 0//min(J_START, J_TARGET)
#define LOWB 0//min(J_START, J_TARGET)
#endif

//Debugging Help
#define CHECK true
#define PRINT true
#define PRINTBEST true

//The change_array variables. 0 -> This method will not be used, #>0 -> use this method on every #th iteration of change_array
#define ROW 1
#define COL 0
#define ALTROW 0
#define ALTCOL 0
#define ALTROW2 0
#define ALTCOL2 0
#define ROWK 0
#define ROWJ 0
#define ROWB 0
#define COLK 0
#define COLJ 0
#define COLB 0
#define SINGLE 0
#define VARIATIONS (ROW && 1)+(COL && 1)+(ALTROW && 1)+(ALTCOL && 1)+(ALTROW2 && 1)+(ALTCOL2 &&  1)+(SINGLE &&  1)+(ROWK && 1)+(ROWJ && 1)+(ROWB && 1)+(COLK && 1)+(COLJ && 1)+(COLB && 1)




using namespace std;


double monte_carlo(int *table, unsigned long long int *b, int num_electrons, int N, int lattice[][NX], double*ham_mod, double* psi_start, double temp, int total_steps, double time_step, gsl_rng * r, double *k_best,double *j_best, double *b_best);
void evolve(int *table, unsigned long long int *b,int num_electrons,int N,  int lattice[][NX], double *psi, double *k_array, double *j_array, double *b_array, int *bonds, int total_steps, double time_step);
double cost(double *psi, double *ham_mod, int N);
double calc_initial_temp(int *table, unsigned long long int *b, int num_electrons, int N, int lattice[][NX], double*ham_mod, int total_steps, double time_step, double *psi_start, gsl_rng *r);
void construct_device_hamiltonian(int *table, unsigned long long int *b,int num_electrons,int N, double *ham_dev, int lattice[][NX],double *k_array, double *j_array, double *b_array, int *bonds, int index, int D);
void construct_model_hamiltonian(int *table, unsigned long long int *b,int num_electrons,int N, double *ham_mod,  int lattice[][NX], int D);
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
void assign_bonds(int *bonds, int lattice[][NX]);
void copy_arrays(int N, double *k_array, double *j_array, double* b_array,  double* k_best,  double* j_best, double* b_best, int total_steps);
void init_arrays(double *k_array, double *j_array,double *b_array, gsl_rng *r, int total_steps);
void change_array(double *k_array, double *j_array, double *b_array, int random_row, int random_col, double change_pm, int i, int total_steps);
void change_row(double *k_array, double *j_array,double *b_array, int row, double change, bool k, bool j, bool b, int jump, int offset);
void change_col(int total_steps,double *k_array,double *j_array,double *b_array, int col, double change,bool k, bool j, bool b, int jump, int offset);
void change_single(double *k_array,double *j_array,double *b_array, int row,int col, double change);
void construct_lattice(int lattice[][NX]);
int get_neighbors(int site, int *neighbors, int lattice[][NX], int D);
double get_random(double lower, double upper, gsl_rng *r);
void print_best_arrays(double *k_array, double *j_array, double *b_array, int totals_step);
void check_norm(double* psi, int N);
void check_unitary(double* ham, int N);
void check_hermicity(double* ham_mod, int N);
void check_weights(double* state, double* ham, int N);
void print_vector(double* psi,int N);
void print_hamiltonian(double* hamiltonian, int N);
void print_hamiltonianR(double *hamiltonian, int N);
void print_E(double ground_E, double best_E);
void test_function();
extern "C" int zgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K, double *ALPHA, double *Z, int *LDA, double *X, int *LDB, double *BETA, double *Y, int *LDC); //complex matrix*matrix mult, odd indices hold imaginary values.
extern "C" int zgemv_(char *TRANS, int *M, int *N,double *ALPHA,double *A, int *LDA, double *X, int *INCX, double *BETA, double *Y, int *INCY); //complex matrix-vector mult, odd indices hold imaginary values.
extern "C" int dsyev_(char *JOBZ, char *UPLO, int *N, double *Vdag, int *LDA, double *D, double *WORK, int *LWORK, int *INFO);//diagonalization, returns the eigenvectors in Vdag and eigenvalues in D.
extern "C" int zgetrf_ (int *M, int *N, double *D, int *LDA, int *IPIV, int *INFO);//A matrix factorization?
extern "C" int zgetrs_( char *TRANS, int *N, int *NRHS,double *D, int *LDA, int *IPIV,double *E, int *LDB, int *INFO);//solves a system of equations
extern "C" void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
extern "C" void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
extern "C" int dgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K, double *ALPHA, double *Z, int *LDA, double *X, int *LDB, double *BETA, double *Y, int *LDC); //real values matrix*matrix mult
extern "C" double zdotc_(int *N, double*ZX,int *INCX, double *ZY, int *INCY);//dots the complex conjugate of ZX with ZY




int main (int argc, char *argv[])
{
	if(PRINT) printf("\nSeed: %i\n", SEED);
	int *table,*bonds,lattice[NY][NX],num_electrons,N,i,j, z, t,x, total_time, total_steps, y, ts, p, multiplier = 2;
	unsigned long long int *b;
	double *ham_target, *ham_start, *k_best, *j_best, *b_best, *k_ground, *j_ground, *b_ground, *psi_start, ground_E, initial_temp=1, best_E, b0,j0,k0, seed = SEED, time_step;
	ofstream file;

	int total_times[] = {2};//{1,2,4};
	int total_steps_array[] = {1,2,4,8,16};


	int max_steps = multiplier*total_steps_array[sizeof(total_steps_array)/sizeof(total_steps_array[0]) - 1]; //getting the last element


	gsl_rng_env_setup();
	const gsl_rng_type * TT = gsl_rng_default;
	gsl_rng * r  = gsl_rng_alloc (TT);
	gsl_rng_set(r, seed);

	construct_lattice(lattice);

	//change this later to cycle through all electrons
	for(i=0;i<1;i++)
	{
		num_electrons=i+1;
		N=choose(num_electrons);

		b = (unsigned long long int*) malloc(N*sizeof(double));
		table=(int*) malloc(num_electrons*N*sizeof(int));
		ham_target = (double *) malloc (2*N*N*sizeof (double));
		psi_start = (double *) malloc(2*N*sizeof(double));
		k_ground = (double *) malloc(NUM_SITES*4*sizeof(double));
		j_ground = (double *) malloc(NUM_SITES*4*sizeof(double));
		b_ground = (double *) malloc(NUM_SITES*2*sizeof(double));
		bonds = (int*) malloc(NUM_SITES*NUM_SITES*sizeof(int));
		ham_start = (double *) malloc(2*N*N*sizeof(double));

		combinations (num_electrons,b,table, N);
		for(j=0;j<NUM_SITES*NUM_SITES;j++) bonds[j]=1;


		//change bounds later to cycle through multiple hams
		for(b0=0.5;b0<0.6;b0+=0.1)
		{
			for(k0=0.5;k0<0.6;k0+=0.1)
			{
				for(j0=0.5;j0<0.6;j0+=0.1)
				{
					j0 += 0.2, k0+= 0.2; b0+=0.3;
					k_ground[0] = K_START, j_ground[0] = J_START, k_ground[NUM_SITES*2] = k0, j_ground[NUM_SITES*2] = j0;
					for(z=0;z<NUM_SITES;z++) b_ground[z] = B_START, b_ground[NUM_SITES+z] = b0;


					construct_device_hamiltonian(table, b, num_electrons, N, ham_start,  lattice, k_ground, j_ground, b_ground, bonds, 0, DEVICE_DIMENSION);
					construct_device_hamiltonian(table, b, num_electrons, N, ham_target, lattice, k_ground, j_ground, b_ground, bonds, 1, DEVICE_DIMENSION);


					get_ground_state(N, ham_start,psi_start);
					ground_E = get_ground_E(N, ham_target);
					if(CHECK) check_norm(psi_start, N);



					k_best = (double *) malloc(2*NUM_SITES*max_steps*sizeof(double));
					j_best = (double *) malloc(2*NUM_SITES*max_steps*sizeof(double));
					b_best = (double *) malloc(NUM_SITES*max_steps*sizeof(double));
					for (x=0;x<2*NUM_SITES*max_steps;x++) k_best[x] =0, j_best[x] = 0;
					for (x=0;x<NUM_SITES*max_steps;x++) b_best[x] =0;
					init_arrays(k_best, j_best, b_best, r, 1);

					for (t = 0; t<sizeof(total_times)/sizeof(total_times[0]);t++)
					{
						total_time =total_times[t];
						file.open("E_VS_T_TOTTIME=" + to_string(total_time) + ".txt");
						for(ts=0;ts<sizeof(total_steps_array)/sizeof(total_steps_array[0]);ts++)
						{
							total_steps = total_steps_array[ts];
							time_step = total_time/((double) total_steps);

							char output[200];
							sprintf(output, "Total Time, Time Step, Total Steps\n %i %f %i\n", total_time, time_step,total_steps);
							file << output;


							initial_temp = calc_initial_temp(table, b, num_electrons, N, lattice, ham_target,total_steps, time_step, psi_start, r);
							best_E = monte_carlo(table, b, num_electrons, N, lattice, ham_target, psi_start, initial_temp,total_steps, time_step,r, k_best, j_best, b_best);


							assign_bonds(bonds, lattice);
							double *E_list, *T_list, *psi, *ham_evolve;
							E_list = (double*) malloc(total_steps*sizeof(double));
							T_list = (double*) malloc(total_steps*sizeof(double));

							ham_evolve = (double *) malloc(2*N*N*sizeof(double));

							psi = (double*) malloc(N*2*sizeof(double));

							memcpy(psi,psi_start, 2*N*sizeof(double));
							for(x=0;x<total_steps;x++)
							{
								construct_device_hamiltonian(table, b, num_electrons, N, ham_evolve, lattice, k_best, j_best, b_best, bonds,x,DEVICE_DIMENSION);
								double *ham_t_i, *exp_matrix;
								exp_matrix = (double *) malloc (2*N*N*sizeof (double));
								ham_t_i = (double *) malloc (2*N*N*sizeof (double));

								for (y=0; y<N*N*2; y++) exp_matrix[y] = 0.0, ham_t_i[y]=0.0;
								for (y=0; y<N*N; y++)
								{
									ham_t_i[2*y+1] = (ham_evolve[2*y]*-time_step); //multiplying by -i*dt for the Pade approximation
									ham_t_i[2*y] = (ham_evolve[2*y+1]*time_step);
								}
								exp_general_complex_double(N, ham_t_i, exp_matrix);
								if(CHECK) check_unitary(exp_matrix, N);
								matrix_vector_mult(exp_matrix, psi, N);

								//print_vector(psi, N);

								if(CHECK) check_norm(psi, N);
								E_list[x] = cost(psi, ham_target, N);
								T_list[x] = (x+1)*time_step;

								free(ham_t_i), free(exp_matrix);
							}

							file << "times:\n[";
							for(x=0;x<total_steps;x++)
							{
								file << ", " << T_list[x];
							}
							file << "]\nExpectation-E_ground:\n[";
							for(x=0;x<total_steps;x++)
							{
								file << ", " << (E_list[x]-ground_E);
							}
							double *k_avg, *b_avg, *j_avg, j_av,k_av,b_av;
							k_avg = (double*) malloc(total_steps*sizeof(double));
							j_avg = (double*) malloc(total_steps*sizeof(double));
							b_avg = (double*) malloc(total_steps*sizeof(double));
							if (PRINTBEST) printf("\nPrinting the best K-J-B arrays from main"), print_best_arrays(k_best, j_best, b_best, total_steps);
							for(x=0;x<total_steps;x++)
							{
								k_av =0, b_av =0, j_av=0;
								for(p=0;p<NUM_SITES*2;p++)
								{
									j_av += j_best[x*NUM_SITES*2+p];
									k_av += k_best[x*NUM_SITES*2+p];
								}
								for(p=0;p<NUM_SITES;p++)
							{
								b_av += b_best[x*NUM_SITES+p];
							}
							k_avg[x] = k_av/(NUM_SITES*2);
							j_avg[x] = j_av/(NUM_SITES*2);
							b_avg[x] = b_av/(NUM_SITES);
							}
							file << "]\nK averages\n[";
							for(x=0;x<total_steps;x++)
							{
								file << ", " << k_avg[x];
							}

							file << "]\nJ averages\n[";
							for(x=0;x<total_steps;x++)
							{
								file << ", " << j_avg[x];
							}
							file << "]\nB averages\n[";
							for(x=0;x<total_steps;x++)
							{
								file << ", " << b_avg[x];
							}
							file << "]\n\n\n";
							if(PRINT) print_E(ground_E, best_E);



//Scaling the arrays by two, using 1 array. pushing the data back
							for(x=total_steps;x>0;x--)
							{
								for(p=0;p<NUM_SITES*2;p++){
									k_best[NUM_SITES*2*(multiplier*x-1)+p] = 	k_best[NUM_SITES*2*(x-1)+p];
								  k_best[NUM_SITES*2*(multiplier*x-2)+p] = 	k_best[NUM_SITES*2*(x-1)+p];
									j_best[NUM_SITES*2*(multiplier*x-1)+p] = 	j_best[NUM_SITES*2*(x-1)+p];
								  j_best[NUM_SITES*2*(multiplier*x-2)+p] = 	j_best[NUM_SITES*2*(x-1)+p];
								}
								for(p=0;p<NUM_SITES;p++){
									b_best[NUM_SITES*(multiplier*x-1)+p] = 	b_best[NUM_SITES*(x-1)+p];
								  b_best[NUM_SITES*(multiplier*x-2)+p] = 	b_best[NUM_SITES*(x-1)+p];
								}
							}
						}
						file.close();
					}
				}
			}
		}
		free(ham_start), free(bonds), free(k_ground), free(j_ground), free(b_ground), free(ham_target), free(b), free(table), free(k_best), free(j_best), free(b_best);
	}
	exit (0);
}



double monte_carlo(int *table, unsigned long long int *b, int num_electrons, int N, int lattice[][NX], double*ham_mod, double* psi_start, double temperature, int total_steps, double time_step, gsl_rng * r,double *k_best,double *j_best, double *b_best)
/*monte_carlo the values in the j, b, and k list in order to produce the lowest energy (expectation value) between the final state (psi), produced by evolve, and the model hamiltonian. This is done by randomly selecting one row of each list, making a slight change, then determining if the new energy is lower than the old. If the new energy is greater than the old, keep with probability exp(delta_E/Temp)*/
{
	int i=0,j=0, random_row=0, random_col, proposal_accepted=0,proposal_count=0, poor_acceptance_count=0,*bonds;
	double *psi, *k_array, *j_array,*b_array,*k_temp, *j_temp, *b_temp, acceptance_rate=0, E_old=0, E_new=0,E_best=0, change_pm=0;
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
	//copy_arrays(N, k_array, j_array, b_array, k_best, j_best,b_best, total_steps);
	evolve(table,b,num_electrons,N,lattice, psi, k_array,j_array,b_array,bonds, total_steps, time_step);
	if(PRINTBEST) print_vector(psi,N);
	if(PRINTBEST) print_hamiltonian(ham_mod, N);
	E_best = cost(psi, ham_mod, N);
	E_old = E_best;
	if(PRINT) printf("\nPre-Monte_Carlo Expectation:   %f\n", E_best);

	for (i=0;i<TEMP_DECAY_ITERATIONS;i++)
	{
		proposal_accepted = 0;
		proposal_count = 0;

		for (j=0; j<SWEEPS;j++)
		{
			copy_arrays(N, k_array, j_array, b_array, k_temp, j_temp,b_temp,total_steps);//a temporary array, used in the undoing of the changes

			change_pm = pow(-1,(int)floor(get_random(0,10,r))) * CHANGE;
			random_row = floor(get_random(0,total_steps,r));
			random_col = floor(get_random(0,NUM_SITES*2,r));

			change_array(k_array,j_array,b_array,random_row,random_col,change_pm,j, total_steps);
			memcpy(psi,psi_start, 2*N*sizeof(double));//resetting psi

			evolve(table,b,num_electrons,N,lattice, psi, k_array,j_array,b_array, bonds, total_steps, time_step);
			E_new = cost(psi, ham_mod, N);
			if (E_new<E_best) E_best=E_new, E_old=E_new,  copy_arrays(N, k_array, j_array, b_array, k_best, j_best,b_best, total_steps);
			else if (E_new<E_old) E_old=E_new;
			else if (get_random(0,1,r)<exp(-(E_new-E_old)/(temperature))) E_old=E_new, proposal_accepted++, proposal_count++;
			else copy_arrays(N, k_temp, j_temp, b_temp, k_array, j_array, b_array,total_steps),proposal_count++;//undoing the change
		}

		acceptance_rate = (double)proposal_accepted/proposal_count;

		if(PRINT) printf("accepted_props:%3i |total_props:%3i |AcceptanceRate: %3.4f |New Expectation: %3.6f\n", proposal_accepted, proposal_count,acceptance_rate,E_old);

		if(acceptance_rate<0.011) poor_acceptance_count++;
		else poor_acceptance_count = 0;

		if(poor_acceptance_count>LIMIT)
		{
			printf("NO PROGRESS FOR %i TEMP CHANGE ITERATIONS, MOVING TO NEXT START STATE\n", LIMIT);
		       	goto end_loop;
		}
		temperature=temperature*EXP_DECAY;
	}
	end_loop:
	if(PRINT) printf("Post-Monte_Carlo Expectation:  %f\n", E_best);
	//if(PRINTBEST) printf("\nBest JKB array values:"), print_best_arrays(k_best, j_best, b_best, total_steps);//print_best_arrays(k_array, j_array, b_array);

	/*memcpy(psi,psi_start,N*2*sizeof(double));
	evolve(table,b,num_electrons,N,lattice, psi, k_best,j_best,b_best,bonds, total_steps, time_step);
	printf("\n\n\nTESTCOST1 = %f\n",cost(psi, ham_mod, N));
*/

	free(k_array), free(j_array), free(b_array),free(k_temp), free(j_temp), free(b_temp),free(psi), free(bonds);
	return E_best;
}




void evolve(int *table, unsigned long long int *b,int num_electrons,int N, int lattice[][NX], double *psi, double *k_array, double *j_array, double *b_array, int *bonds, int total_steps, double time_step)
/*evolve a starting state, psi, by acting on it with exp(ham_dev*-i*time_step). The resulting state is updated as psi and the process repeats total_steps times (total_time/total_steps) until the final psi state is produced. The function contains two methods for calculating exp(ham_dev*-i+time_step), one is a diagonalization method, the other a Pade approximation*/
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
	}
	free(ham_dev), free(ham_t_i), free(exp_matrix), free(D), free(Vdag),free(ham_real);
}




double cost(double *psi, double *ham_mod, int N)
/*Computing the expectation value between ham_mod and psi, <psi|ham_mod|psi>*/
{

	int i=0,INCX = 1,INCY = 1,j;
	double *psi_conj, *psi_temp, result=0, resulti=0, norm;
	psi_conj = (double*) malloc (N*2*sizeof(double));
	psi_temp = (double*) malloc (N*2*sizeof(double));
	memcpy(psi_conj, psi, 2*N*sizeof(double));
	memcpy(psi_temp, psi, 2*N*sizeof(double));

	if(CHECK) check_norm(psi_temp, N);
	if(CHECK) check_weights(psi_temp, ham_mod, N);

	matrix_vector_mult(ham_mod, psi_temp, N);//H*psi, the operator acting on the ket, storing in psi
	result = zdotc_(&N, psi_conj, &INCX, psi_temp, &INCY);//psi* * psi, the dot between the complex conj and the result of H*psi

	free(psi_conj), free(psi_temp);
//	print_hamiltonian(ham_mod, N);
	//print_vector(psi, N);

	//if(result > -1.59) printf("%f\n", result);
	return result;
}




double calc_initial_temp(int *table, unsigned long long int *b, int num_electrons, int N, int lattice[][NX], double*ham_target, int total_steps, double time_step, double *psi_start, gsl_rng *r)
/*Finding an good initial temp that allows an average increase acceptance probability of about ACCEPTANCE_PROB (typically 0.8). Choosing an initial temperature that allows the chance of accepting j,b, and k_array values which increase the expectation value to be ~80%, https://www.phy.ornl.gov/csep/mo/node32.html*/
{
	if(PRINT) printf("\n\n\n\n\n\n\n\n\n\n\n\n\n...Calculating initial temperature based on %i random starting states...\n", RANDOM_STATES);
	int *bonds, i=0,j=0, random_row=0,random_col=0, start_state=0, count=0;
	double *psi,*psi_random, *k_array, *j_array,*b_array, E_old=0, E_new=0,sum=0, change_pm=0, initial_temp=0;
	bonds = (int*) malloc(NUM_SITES*NUM_SITES*sizeof(int));
	k_array = (double *) malloc(2*NUM_SITES*total_steps*sizeof(double));
	j_array = (double *) malloc(2*NUM_SITES*total_steps*sizeof(double));
	b_array = (double *) malloc(NUM_SITES*total_steps*sizeof(double));
	psi_random = (double *) malloc (2*N*sizeof(double));
	psi = (double *) malloc (2*N*sizeof(double));

	assign_bonds(bonds, lattice);
	for (i=0; i<N*2;i++) psi[i]=0.0;
	for (j=0;j<RANDOM_STATES;j++)
	{
		for (i=0; i<N*2;i++) psi_random[i] =0.0;
		start_state = floor(get_random(0,N,r));

		psi_random[start_state*2] = 1;
		memcpy(psi,psi_random, 2*N*sizeof(double));

		(k_array, j_array, b_array,r, total_steps);

		evolve(table,b,num_electrons,N,lattice, psi, k_array,j_array,b_array, bonds, total_steps, time_step);
		E_old = cost(psi, ham_target, N);

		for (i=0; i<SWEEPS;i++)
		{
			change_pm = pow(-1,(int)floor(get_random(0,10,r))) * CHANGE;
			random_row = floor(get_random(0,total_steps,r));
			random_col = floor(get_random(0,NUM_SITES*2,r));

			change_array(k_array,j_array,b_array,random_row,random_col,change_pm,i, total_steps);
			memcpy(psi,psi_random, 2*N*sizeof(double));//resetting psi
			evolve(table,b,num_electrons,N,lattice, psi, k_array,j_array,b_array, bonds, total_steps, time_step);
			E_new = cost(psi, ham_target, N);

			if (E_new>E_old) sum += (E_new-E_old), count++;
			E_old=E_new;
		}
	}
	free(k_array), free(j_array), free(b_array), free(psi), free(psi_random), free(bonds);
	initial_temp = -(sum/(count*log(ACCEPTANCE_PROB)));

	if(PRINT){
		printf("#######################################################################");
		printf("\nELECTRONS:   %i\nDIMENSION:   %i\nTOTAL_TIME:  %i\nTOTAL_STEPS: %i\nTIME_STEP:   %f\nINIT_TEMP:   %f\n",num_electrons, N,(int) time_step*total_steps,total_steps, time_step,initial_temp);
		print_vector(psi_start, N);
		printf("#######################################################################\n");
	}

	return initial_temp;
}




void construct_device_hamiltonian(int *table, unsigned long long int *b,int num_electrons,int N, double *ham_dev, int lattice[][NX],double *k_array, double *j_array, double *b_array, int *bonds,  int index, int D)
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
	free(neighbors);
	if(CHECK) check_hermicity(ham_dev, N);
}




void construct_model_hamiltonian(int *table, unsigned long long int *b,int num_electrons,int N, double *ham_mod, int lattice[][NX], int D)
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
	free(neighbors);
	if(CHECK) check_hermicity(ham_mod, N);
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
	if (INFO !=0) printf("DIAGONALIZATION ERROR, INFO = %i\n", INFO);
       	LWORK=WORK[0];
       	free(WORK);
       	WORK=(double*) malloc(LWORK*sizeof(double));

       	dsyev_(&JOBZ, &UPLO, &N, Vdag, &LDA, D, WORK, &LWORK, &INFO );
       	if (INFO !=0) printf("DIAGONALIZATION ERROR, INFO = %i\n", INFO);
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
	if (INFO !=0) printf("ERROR, INFO = %i\n", INFO);
	zgetrs_( &TRANS, &N, &NRHS,D,    &LDA, IPIV,E, &LDB, &INFO);//solving a system of equations
	if (INFO !=0) printf("ERROR, INFO = %i\n", INFO);

	for (kk=1;kk<=s;kk++)
	{
		memcpy(X,E,2*N*N*sizeof(double));
		zgemm_ (&TRANSA, &TRANSB, &M, &N, &K, ALPHA, E, &LDA, X, &LDB, BETA, Y, &LDC);//matrixmultiplication
		if (INFO !=0) printf("ERROR, INFO = %i\n", INFO);
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
	for(i=0; i<N; i++) if(D[i]<min_E){ min_E = D[i]; min_i = i;}
	for(i=0; i<N; i++) ground_state[i*2] = Vdag[min_i*N+i];

	free(Vdag), free(D), free(ham_real);
}




double get_ground_E(int N, double *ham)
/*Finding the ground energy, the smallest eigenvalue of the matrix*/
{
	print_hamiltonian(ham, N);
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




void assign_bonds(int *bonds, int lattice[][NX])
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
}




void copy_arrays(int N, double *k_array, double *j_array, double* b_array,  double* k_to,  double* j_to, double* b_to, int total_steps)
/*storing the updated k, j, and b values*/
{
	memcpy(k_to, k_array, 2*NUM_SITES*total_steps*sizeof(double));
	memcpy(j_to, j_array, 2*NUM_SITES*total_steps*sizeof(double));
	memcpy(b_to, b_array, NUM_SITES*total_steps*sizeof(double));
}




void init_arrays(double *k_array,double *j_array,double *b_array,gsl_rng * r, int total_steps)
/*Initializng the values of the k, j, and b lists which hold the values of the constants for each site (b_array) and between each bond (k_array and j_array)*/
{
	int i,j;
	double upj=1.0,lowj=-1.0,upk=1.0,lowk=-1.0,upb=1.0,lowb=-1.0, random_val1, random_val2, random_val3;
	if(DEVICE_SAME_MODEL) upj=UPJ,lowj=LOWJ,upk=UPK,lowk=LOWK,upb=UPB,lowb=LOWB;
	if(UNIFORM)
	{
		for (i=0; i<total_steps;i++)
		{

			random_val1 = get_random(lowj,upj,r);
			random_val2 = get_random(lowk,upk,r);
			random_val3 = get_random(lowb,upb,r);
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
				j_array[i*NUM_SITES*2+j] = get_random(lowj,upj,r);
				k_array[i*NUM_SITES*2+j] = get_random(lowk,upk,r);
			}
			for (j=0; j<NUM_SITES; j++)
			{
				b_array[i*NUM_SITES+j]= get_random(lowb,upb,r);
			}
		}
	}
	if(PRINTBEST) printf("\nPrinting the best K-J-B arrays from init_arrays"), print_best_arrays(k_array, j_array, b_array, total_steps);
}




void change_array(double *k_array, double *j_array, double *b_array, int random_row, int random_col, double change_pm, int i, int total_steps)
/*changing the j, k, and b arrays. Use the defines at start of program to determine which manipulation functions will be used, and in what order*/
{
	int mod = VARIATIONS;
	if(i%mod==ROW-1) change_row(k_array,j_array,b_array,random_row,change_pm,true, true,true, 1, 0);
	else if(i%mod==COL-1) change_col(total_steps,k_array,j_array,b_array,random_col,change_pm, true, true, true, 1, 0);
	else if(i%mod==ALTROW-1) change_row(k_array,j_array,b_array,random_row,change_pm, true, true, true ,2, 0);
	else if(i%mod==ALTCOL-1) change_col(total_steps,k_array,j_array,b_array,random_col,change_pm, true, true, true, 2, 0);
	else if(i%mod==ALTROW2-1) change_row(k_array,j_array,b_array,random_row,change_pm, true, true, true, 2, 1);
	else if(i%mod==ALTCOL2-1) change_col(total_steps,k_array,j_array,b_array,random_col,change_pm, true, true, true, 2, 1);
	else if(i%mod==ROWK-1) change_row(k_array,j_array,b_array,random_row,change_pm, true, false, false, 1, 0);
	else if(i%mod==ROWJ-1) change_row(k_array,j_array,b_array,random_row,change_pm, false, true, false, 1, 0);
	else if(i%mod==ROWB-1) change_row(k_array,j_array,b_array,random_row,change_pm, false, false, true, 1, 0);
	else if(i%mod==COLK-1) change_col(total_steps,k_array,j_array,b_array,random_col,change_pm, true, false, false, 2, 1);
	else if(i%mod==COLJ-1) change_col(total_steps,k_array,j_array,b_array,random_col,change_pm, false, true, false, 2, 1);
	else if(i%mod==COLB-1) change_col(total_steps,k_array,j_array,b_array,random_col,change_pm, false, false, true, 2, 1);
	else if(i%mod==SINGLE-1) change_single(k_array,j_array,b_array,random_row,random_col,change_pm);
	else printf("NOPE\n");
}



void change_row(double *k_array,double *j_array,double *b_array, int row, double change, bool k, bool j, bool b, int jump, int offset)
/*Changing all of the lists by a value change at at the row number, row. Used in the monte_carlo function. offset gives the starting element, jump gives the amount to increase each increment.
  bool k, j, and b determine if the k,j,and b lists will be changes. Bounds of the values that the lists can take are given by the #defines if using a device that's the same as model, +-1 otherwise*/
{
	int i;
	double upj=1.0,lowj=-1.0,upk=1.0,lowk=-1.0,upb=1.0,lowb=-1.0;
	if(DEVICE_SAME_MODEL) upj=UPJ,lowj=LOWJ,upk=UPK,lowk=LOWK,upb=UPB,lowb=LOWB;


  if(k) for (i=offset; i<NUM_SITES*2; i+=jump)
	{
		//printf("ROW: %i", row);
		//printf("NUM_SITES*2*row+i: %i\n", NUM_SITES*2*row+i);
		if ((lowk < k_array[NUM_SITES*2*row+i] + change) &&  (k_array[NUM_SITES*2*row+i] + change < upk)) k_array[NUM_SITES*2*row+i] += change;
	}
	if(j) for (i=offset; i<NUM_SITES*2; i+=jump)
	{
		if ((lowj < j_array[NUM_SITES*2*row+i] + change) && (j_array[NUM_SITES*2*row+i] + change < upj)) j_array[NUM_SITES*2*row+i] += change;
	}
	if(b) for (i=offset; i<NUM_SITES; i+=jump)
	{
		if ((lowb < b_array[NUM_SITES*row+i] + change) && (b_array[NUM_SITES*row+i] + change < upb)) b_array[NUM_SITES*row+i] += change;
	}

	//if(PRINTBEST) printf("Best JKB array values:"), print_best_arrays(k_array, j_array, b_array,16 );
}




void change_col(int total_steps,double *k_array,double *j_array,double *b_array, int col, double change,bool k, bool j, bool b, int jump, int offset)
/*Changing all of the lists by a value change at at the col number, col. Used in the monte_carlo function. offset gives the starting element, jump gives the amount to increase each increment.
  bool k, j, and b determine if the k,j,and b lists will be changes. Bounds of the values that the lists can take are given by the #defines if using a device that's the same as model, +-1 otherwise*/
{
	int i;
	double upj=1.0,lowj=-1.0,upk=1.0,lowk=-1.0,upb=1.0,lowb=-1.0;
	if(DEVICE_SAME_MODEL) upj=UPJ,lowj=LOWJ,upk=UPK,lowk=LOWK,upb=UPB,lowb=LOWB;

	if(k) for (i=offset; i<total_steps; i+=jump)
	{
		if ((lowk < k_array[NUM_SITES*2*i+col] + change) && (k_array[NUM_SITES*2*i+col] + change  < upk)) k_array[NUM_SITES*2*i+col] += change;
	}
	if(j) for (i=offset; i<total_steps; i+=jump)
	{
		if ((lowj < j_array[NUM_SITES*2*i+col] + change) && (j_array[NUM_SITES*2*i+col] + change < upj)) j_array[NUM_SITES*2*i+col] += change;
	}
	if(b) for (i=offset; i<total_steps; i+=jump)
	{
		if ((lowb < b_array[NUM_SITES*i+(int)floor(col/2.0)] + change) && (b_array[NUM_SITES*i+(int)floor(col/2.0)] + change < upb)) b_array[NUM_SITES*i+(int)floor(col/2.0)] += change;
	}
}




void change_single(double *k_array,double *j_array,double *b_array, int row,int col, double change)
/*Changing a single element in each array, at column col, row row*/
{
	double upj=1.0,lowj=-1.0,upk=1.0,lowk=-1.0,upb=1.0,lowb=-1.0;
	if(DEVICE_SAME_MODEL) upj=UPJ,lowj=LOWJ,upk=UPK,lowk=LOWK,upb=UPB,lowb=LOWB;

	if ((lowj < j_array[NUM_SITES*2*row+col] + change) && (j_array[NUM_SITES*2*row+col] + change < upj)) j_array[NUM_SITES*2*row+col] += change;
 	if ((lowk < k_array[NUM_SITES*2*row+col] + change) && (k_array[NUM_SITES*2*row+col] + change < upk)) k_array[NUM_SITES*2*row+col] += change;
	if ((lowb < b_array[NUM_SITES*row+(int)floor(col/2.0)] + change) && (b_array[NUM_SITES*row+(int)floor(col/2.0)] + change  < upb)) b_array[NUM_SITES*row+(int)floor(col/2.0)] += change;
}




void construct_lattice(int lattice[][NX])
/*Build a NX*NX lattice, with sites 1 through NX*NX listed in a snaking pattern*/
{
	int x,y;
	for (x=0;x<NX;x++)
	{
		if (x%2 ==1) for (y=0;y<NX;y++) lattice[x][y] = (NX*x)+y+1;
		else for (y=0;y<NX;y++) lattice[x][NX-1-y] = NX*x+y+1;
	}
}




int get_neighbors(int site, int *neighbors, int lattice[][NX], int D)
/*Gets the neighbors of a site, returning all 4 neighbors if open_boundry coditions are true, otherwise just closed-boudnary neighbors*/
{
	int x,y,i,count=0;
	if (D==2)
	{
		for (x=0;x<NX;x++) for (y=0;y<NX;y++) if(site == lattice[x][y]) goto end_loop;//Finding the coordinates
		end_loop: for(i=0;i<4;i++) neighbors[i]=0;

		if (!OPEN)
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
			if (y+1<NX) neighbors[count] = lattice[x][y+1], count++;
			if (y>0) neighbors[count] = lattice[x][y-1], count++;
		}
		return count;
	}
	else if (D==1)
	{
		if (OPEN)
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
	printf("ERROR! NO NEIGHBORS FOUND.");
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





void print_best_arrays(double *k_array, double *j_array, double *b_array, int total_steps)
{
	int i,j;

	printf("\nK_array: %iX%i (stepsXsites)\n[", total_steps, NUM_SITES*2);
	for(i=0;i<total_steps;i++)
	{
		for(j=0;j<2*NUM_SITES;j++)
		{
			printf("%5.2f ",k_array[i*NUM_SITES*2+j]);
			if(k_array[i*NUM_SITES+j] > 1 || -1 > k_array[i*NUM_SITES+j]) printf("\n\n\n ERROR IN PRINT BEST \n\n\n");
		}
		if (i==total_steps-1) printf("]\n");
		else printf(";\n ");
	}

	printf("J_array: %iX%i (stepsXsites)\n[", total_steps, NUM_SITES*2);
	for(i=0;i<total_steps;i++)
	{
		for(j=0;j<2*NUM_SITES;j++)
		{
			printf("%5.2f ",j_array[i*NUM_SITES*2+j]);
			if(k_array[i*NUM_SITES+j] > 1 || -1 > k_array[i*NUM_SITES+j]) printf("\n\n\n ERROR IN PRINT BEST  \n\n\n");
		}
		if (i==total_steps-1) printf("]\n");
		else printf(";\n ");
	}

	printf("B_array: %iX%i (stepsXsites)\n[", total_steps, NUM_SITES);
	for(i=0;i<total_steps;i++)
	{
		for(j=0;j<NUM_SITES;j++)
		{
			printf("%5.2f ",b_array[i*NUM_SITES+j]);
			if(b_array[i*NUM_SITES+j] > 1 || -1 > b_array[i*NUM_SITES+j]) printf("\n\n\n ERROR IN PRINT BEST \n\n\n");
		}
		if (i==total_steps-1) printf("]\n");
		else printf(";\n ");
	}
	printf("\n");
}




void check_norm(double* psi, int N)
{
	int i,j=0;
	double sum=0;
	for(i=0;i<N*2; i+=2)
	{
		sum+= psi[i]*psi[i]+(psi[i+1]*psi[i+1]);

	}
	if(sum>1.000000001 or sum<0.99999999) printf("NORM ERROR, SIZE: %f\n", sum);
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
		if(unitary[i] < -0.00000001 or unitary[i] > 0.00000001) printf("ERROR, NON UNITARY ELEMENTS AT %i, VALUE: %f\n", i, unitary[i]);
	}
}





void check_hermicity(double* ham, int N)
{
	int i,j;
	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{
			if(ham[2*(j*N+i)] != ham[2*(i*N+j)]) printf("ERROR, NON HERMITIAN ELEMENT AT i,j: %i,%i\nElementi = %10.7f\n j = %10.7f\n", i*N+j, j*N+i, ham[i*N+j], ham[j*N+i]);
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
	if(c_squared_sum>1.00001 or c_squared_sum <0.999999) printf("ERROR, BAD WEIGHTS\n");
	free(Vdag), free(D), free(ham_real);
}


void print_vector(double* psi,int N)
{
	int i;
	printf("\nPrinting the %i dimensional vector:\n[", N);
	for(i=0;i<N;i++)
	{
		printf("%09.5f+%09.5fi  ", psi[2*i], psi[2*i+1]);
	}
	printf("]\n");
}




void print_hamiltonianR(double* hamiltonian, int N)
/*prints a real-valued hamiltonian in matrix form*/
{
	int i,j;
	printf("\nPrinting the %ix%i Hamiltonian:\n[", N,N);
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
	printf("\nPrinting the %ix%i Hamiltonian:\n[", N,N);
	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++) printf("%09.5f+%09.5fi  ",(hamiltonian[2*(j*N+i)]+0.0), hamiltonian[2*(j*N+i)+1]);
		if (i==N-1) printf("]");
		else printf(";\n ");
	}
	printf("\n");
}




void print_E(double ground_E, double best_E)
{
	printf("#######################################################################\n");
	printf("TARGET ENERGY:  %9.6f\n",ground_E);
	printf("BEST ENERGY:    %9.6f\n",best_E);
	printf("DIFFERENCE:     %9.6f\n",best_E-ground_E);
	printf("#######################################################################\n");
}




void test_function()
{
	int i,j;
	int N = 3;
	double *D, *Vdag, *ham_mod, min_E=1000;
	Vdag = (double*) malloc(N*N*sizeof(double));
	D = (double*) malloc(N*sizeof(double));
}
/*Saving the block from main

int main (int argc, char *argv[])
{
	if(PRINT) printf("\nSeed: %i\n", SEED);
	int *table,*bonds,lattice[NY][NX],num_electrons,N,i,j,iterations;
	unsigned long long int *b;
	double *ham_target, *ham_start, *k_ground, *j_ground, *b_ground, *psi_start, ground_E = 100, initial_temp=1;

	int dimensions[] = {2,3,4,5};

	construct_lattice(lattice);
	gsl_rng_env_setup();
	double seed = SEED;
	const gsl_rng_type * TT = gsl_rng_default;
	gsl_rng * r  = gsl_rng_alloc (TT);
	gsl_rng_set(r, seed);


	for(i=0;i<NUM_SITES-1;i++)
	{

		num_electrons=i+1;
		N=choose(num_electrons);
		b = (unsigned long long int*) malloc(N*sizeof(double));
		table=(int*) malloc(num_electrons*N*sizeof(int));
		ham_target = (double *) malloc (2*N*N*sizeof (double));
		psi_start = (double *) malloc(2*N*sizeof(double));

		combinations (num_electrons,b,table, N);

		if(DEVICE_SAME_MODEL)
		{
			k_ground = (double *) malloc(NUM_SITES*4*sizeof(double));
			j_ground = (double *) malloc(NUM_SITES*4*sizeof(double));
			b_ground = (double *) malloc(NUM_SITES*2*sizeof(double));
			bonds = (int*) malloc(NUM_SITES*NUM_SITES*sizeof(int));
			ham_start = (double *) malloc(2*N*N*sizeof(double));

			for(j=0;j<NUM_SITES*NUM_SITES;j++) bonds[j]=1;
			k_ground[0] = K_START, j_ground[0] = J_START, k_ground[NUM_SITES*2] = K_TARGET, j_ground[NUM_SITES*2] = J_TARGET;
			for(j=0;j<NUM_SITES;j++) b_ground[j] = B_START, b_ground[NUM_SITES+j] = B_TARGET;

			construct_device_hamiltonian(table, b, num_electrons, N, ham_start, lattice, k_ground, j_ground, b_ground, bonds, 0, DEVICE_DIMENSION);
			construct_device_hamiltonian(table, b, num_electrons, N, ham_target , lattice, k_ground, j_ground, b_ground, bonds, 1, DEVICE_DIMENSION);

			get_ground_state(N, ham_start,psi_start);
			ground_E = get_ground_E(N, ham_target);
			if(CHECK) check_norm(psi_start, N);

			initial_temp = calc_initial_temp(table, b, num_electrons, N, lattice, ham_target,r);
			monte_carlo(table, b, num_electrons, N, lattice, ham_target, psi_start, initial_temp,r);

			if(PRINT) print_E(ground_E);
			free(ham_start), free(bonds), free(k_ground), free(j_ground), free(b_ground);
		}
		else
		{
			int start_states[] = {0,N/3,2*N/3};
			construct_model_hamiltonian(table, b, num_electrons, N, ham_target, lattice, 2);
			ground_E =get_ground_E(N, ham_target);

			initial_temp = calc_initial_temp(table, b, num_electrons, N, lattice, ham_target,r);
			for(j=0;j<sizeof(start_states)/sizeof(start_states[0]);j++)
			{
				for (i=0; i<N*2;i++) psi_start[i] =0.0;
				if(PRINT) printf("\nStart State:                 %i\n", start_states[j]);
				psi_start[start_states[j]] = 1;
				monte_carlo(table, b, num_electrons, N, lattice, ham_target, psi_start, initial_temp,r);
			}
			if(PRINT) print_E(ground_E);
		}
		free(b), free(table), free(ham_target), free(psi_start);
	}
	exit (0);
}

*/
