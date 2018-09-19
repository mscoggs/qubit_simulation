/*Questions:
-expectation overshot ground state, GS=-1.8, cost=-1.82

TODO:
-try a handful of different times
-verify that acceptance prob shouldn't weigh in the -delta_e changes
-double check annealing code
-find how long it take to miniminze
*/


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

#define DEVICE_SAME_MODEL true
#define SEED 1
#define PRINT true
#define PRINTBEST false
#define OPEN true
#define NX 4
#define NY NX
#define NUM_SITES NY*NX
#define DEVICE_DIMENSION 2
#define DIAG false
#define T 1
#define V 2
#define J1 -0.9
#define K1 -0.7
#define B1 -0.8
#define J2 0.7
#define K2 0.9
#define B2 0.9
#define LIMIT 30
#define TIME_STEP .5
#define TOTAL_TIME 10 
#define TOTAL_STEPS TOTAL_TIME/TIME_STEP
#define RANDOM_STATES 3
#define SWEEPS 200
#define CHANGE 0.02
#define ACCEPTANCE_PROB 0.6
#define EXP_DECAY 0.99
#define TEMP_DECAY_ITERATIONS 600
/*The following are the change_array variables. 0 -> This method will not be used, #>0 -> use this method on every #th iteration of change_array*/
#define ROW 1
#define COL 2
#define ALTROW 3
#define ALTCOL 4
#define ALTROW2 5
#define ALTCOL2 6
#define SINGLE 7
#define VARIATIONS (ROW && 1)+ (COL && 1) + (ALTROW && 1) + (ALTCOL && 1) + (ALTROW2 && 1) + (ALTCOL2 &&  1) + (SINGLE &&  1)




using namespace std;


void optimize(int *table, unsigned long long int *b, int num_electrons, int N, int lattice[][NX], double*ham_mod, double* psi_start, double temp, gsl_rng * r);
void evolve(int *table, unsigned long long int *b,int num_electrons,int N,  int lattice[][NX], double *psi, double *k_array, double *j_array, double *b_array, int *bonds);
double cost(double *psi, double *ham_mod, int N);
double calc_initial_temp(int *table, unsigned long long int *b, int num_electrons, int N, int lattice[][NX], double*ham_mod, gsl_rng *r);
void construct_device_hamiltonian(int *table, unsigned long long int *b,int num_electrons,int N, double *ham_dev, int lattice[][NX],double *k_array, double *j_array, double *b_array, int *bonds, int index, int D);
void construct_model_hamiltonian(int *table, unsigned long long int *b,int num_electrons,int N, double *ham_mod,  int lattice[][NX], int D);
void diag_hermitian_real_double(int N,  double *A, double *Vdag,double *D);
void exp_diaganolized_mat(double *ham_real, double *Vdag, double* D, int N);
void exp_general_complex_double(int N, double *A, double *B);
void matrix_vector_mult(double *exp_matrix, double *psi, int N);
int hop(unsigned long long int b, unsigned long long int *v,int n, int j);
int find(int N,unsigned long long int *v, unsigned long long int *b);
void get_ground_state(int N, double *ham, double* ground_state);
double get_ground_E(int N, double *ham);
unsigned long long int choose(int num_electrons);
int combinations ( int num_electrons, unsigned long long int *b,int *tab, int N);
void assign_bonds(int *bonds, int lattice[][NX]);
void copy_arrays(int N, double *k_array, double *j_array, double* b_array,  double* k_best,  double* j_best, double* b_best);
void init_arrays(double *k_array, double *j_array,double *b_array, gsl_rng *r);
void change_array(double *k_array, double *j_array, double *b_array, int random_row, int random_col, double change_pm, int i);
void change_row(double *k_array, double *j_array,double *b_array, int row, double change, int jump, int offset);
void change_col(double *k_array,double *j_array,double *b_array, int col, double change, int jump, int offset);
void change_single(double *k_array,double *j_array,double *b_array, int row,int col, double change);
void construct_lattice(int lattice[][NX]);
int get_neighbors(int site, int *neighbors, int lattice[][NX], int D);
double get_random(double lower, double upper, gsl_rng *r);
void print_best_arrays(double *k_array, double *j_array, double *b_array);
void print_hamiltonian(double* hamiltonian, int N, bool print_imaginary);
void print_hamiltonianR(double *hamiltonian, int N);
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
//	test_function();
	int *table,*bonds,lattice[NY][NX],num_electrons,N,i,j,iterations;
	unsigned long long int *b;
	double *ham_mod, *ham_ground, *k_ground, *j_ground, *b_ground, *psi_start, ground_E = 100, initial_temp=1;
	int electron_counts[] = {1};
	construct_lattice(lattice);

	gsl_rng_env_setup();
	double seed = SEED;
	const gsl_rng_type * TT = gsl_rng_default;
	gsl_rng * r  = gsl_rng_alloc (TT);
	gsl_rng_set(r, seed);

	for(i=0;i<sizeof(electron_counts)/sizeof(int);i++)
	{
		num_electrons=electron_counts[i];
		N=choose(num_electrons);
		b = (unsigned long long int*) malloc(N*sizeof(double));
		table=(int*) malloc(num_electrons*N*sizeof(int));
		ham_mod = (double *) malloc (2*N*N*sizeof (double));
		psi_start = (double *) malloc(2*N*sizeof(double));

		combinations (num_electrons,b,table, N);

		if(DEVICE_SAME_MODEL)
		{
			k_ground = (double *) malloc(NUM_SITES*4*sizeof(double));
			j_ground = (double *) malloc(NUM_SITES*4*sizeof(double));
			b_ground = (double *) malloc(NUM_SITES*2*sizeof(double));
			bonds = (int*) malloc(NUM_SITES*NUM_SITES*sizeof(int));
			ham_ground = (double *) malloc (2*N*N*sizeof (double));

			for(j=0;j<NUM_SITES*NUM_SITES;j++) bonds[j]=1;
			k_ground[0] = K1;
			k_ground[NUM_SITES*2] = K2;
			j_ground[0] = J1;
			j_ground[NUM_SITES*2] = J2;
			for(j=0;j<NUM_SITES;j++) b_ground[j] = B1;
			for(j=0;j<NUM_SITES;j++) b_ground[NUM_SITES+j] = B2;

			construct_device_hamiltonian(table, b, num_electrons, N, ham_ground, lattice, k_ground, j_ground, b_ground, bonds, 0, DEVICE_DIMENSION);
			construct_device_hamiltonian(table, b, num_electrons, N, ham_mod, lattice, k_ground, j_ground, b_ground, bonds, 1, DEVICE_DIMENSION);
			get_ground_state(N, ham_ground,psi_start);
			ground_E = get_ground_E(N, ham_mod);
			printf("\n\n#######################################################################\n");
			printf("DEVICE AS MODEL GROUND STATE ENERGY: %f\n",ground_E);
			printf("#######################################################################\n");
			initial_temp = calc_initial_temp(table, b, num_electrons, N, lattice, ham_mod,r);

			optimize(table, b, num_electrons, N, lattice, ham_mod, psi_start, initial_temp,r);
			printf("#######################################################################\n");
			printf("DEVICE AS MODEL GROUND STATE ENERGY: %f\n",ground_E);
			printf("#######################################################################\n");
			free(ham_ground), free(bonds), free(k_ground), free(j_ground), free(b_ground);
		}
		else
		{
			int start_states[] = {N/3,N/3,N/3,N/3};//,2*N/3,2*N/3,2*N/3,2*N/3};
			iterations = sizeof(start_states)/sizeof(start_states[0]);
			construct_model_hamiltonian(table, b, num_electrons, N, ham_mod, lattice, 2);
			ground_E =get_ground_E(N, ham_mod);
			printf("\n\n#######################################################################\n");
			printf("MODEL GROUND STATE ENERGY: %f\n",ground_E);
			printf("#######################################################################\n");
			initial_temp = calc_initial_temp(table, b, num_electrons, N, lattice, ham_mod,r);

			for(j=0;j<iterations;j++)
			{
				for (i=0; i<N*2;i++) psi_start[i] =0.0;
				printf("\nStart State:                 %i\n", start_states[j]);
				psi_start[start_states[j]] = 1;
				optimize(table, b, num_electrons, N, lattice, ham_mod, psi_start, initial_temp,r);
			}			
			printf("#######################################################################\n");
			printf("MODEL GROUND STATE ENERGY: %f\n",ground_E);
			printf("#######################################################################\n");
		}
		free(b), free(table), free(ham_mod), free(psi_start);
	}
	exit (0);
}




void optimize(int *table, unsigned long long int *b, int num_electrons, int N, int lattice[][NX], double*ham_mod, double* psi_start, double temperature, gsl_rng * r )
/*Optimize the values in the j, b, and k list in order to produce the lowest energy (expectation value) between the final state (psi), produced by evolve, and the model hamiltonian. This is done by randomly selecting one row of each list, making a slight change, then determining if the new energy is lower than the old. If the new energy is greater than the old, keep with probability exp(delta_E/Temp)*/
{
	int i=0,j=0, random_row=0, random_col, change_accepted=0,bad_change=0, poor_acceptance_count=0,*bonds;
	double *psi, *k_array, *j_array,*b_array,*k_best, *j_best, *b_best, *k_temp, *j_temp, *b_temp, acceptance_rate=0, E_old=0, E_new=0,E_best=0, change_pm=0;
	bonds = (int*) malloc(NUM_SITES*NUM_SITES*sizeof(int));
	k_array = (double *) malloc(2*NUM_SITES*TOTAL_STEPS*sizeof(double));
	j_array = (double *) malloc(2*NUM_SITES*TOTAL_STEPS*sizeof(double));
	b_array = (double *) malloc(NUM_SITES*TOTAL_STEPS*sizeof(double));
	k_best = (double *) malloc(2*NUM_SITES*TOTAL_STEPS*sizeof(double));
	j_best = (double *) malloc(2*NUM_SITES*TOTAL_STEPS*sizeof(double));
	b_best = (double *) malloc(NUM_SITES*TOTAL_STEPS*sizeof(double));
	k_temp = (double *) malloc(2*NUM_SITES*TOTAL_STEPS*sizeof(double));
	j_temp = (double *) malloc(2*NUM_SITES*TOTAL_STEPS*sizeof(double));
	b_temp = (double *) malloc(NUM_SITES*TOTAL_STEPS*sizeof(double));
	psi = (double *) malloc (2*N*sizeof(double));

	assign_bonds(bonds, lattice);
	memcpy(psi,psi_start, 2*N*sizeof(double));

	init_arrays(k_array, j_array, b_array, r);
	copy_arrays(N, k_array, j_array, b_array, k_best, j_best,b_best);
	evolve(table,b,num_electrons,N,lattice, psi, k_array,j_array,b_array,bonds);
	E_best = cost(psi, ham_mod, N);
	E_old = E_best;
	printf("Pre-Optimized Expectation:   %f\n", E_best);

	for (i=0;i<TEMP_DECAY_ITERATIONS;i++)
	{
		change_accepted = 0;
		bad_change = 0;

		for (j=0; j<SWEEPS;j++)
		{
			copy_arrays(N, k_array, j_array, b_array, k_temp, j_temp,b_temp);//a temporary array, used in the undoing of the changes
			change_pm = pow(-1,(int)floor(get_random(0,10,r))) * CHANGE;
			random_row = floor(get_random(0,TOTAL_STEPS,r));
			random_col = floor(get_random(0,NUM_SITES*2,r));
			change_array(k_array,j_array,b_array,random_row,random_col,change_pm,j);
			memcpy(psi,psi_start, 2*N*sizeof(double));//resetting psi
			evolve(table,b,num_electrons,N,lattice, psi, k_array,j_array,b_array, bonds);
			E_new = cost(psi, ham_mod, N);
			if (E_new<E_best) E_best=E_new, copy_arrays(N, k_array, j_array, b_array, k_best, j_best,b_best);

			if (E_new<E_old) E_old=E_new;
			else if (get_random(0,1,r)<exp(-(E_new-E_old)/(temperature))) E_old=E_new, change_accepted++, bad_change++;
			else copy_arrays(N, k_temp, j_temp, b_temp, k_array, j_array, b_array),bad_change++;//undoing the change
		}

		acceptance_rate = (double)change_accepted/bad_change;
		if(PRINT) printf("accepted_props:%3i |total_props:%3i |AcceptanceRate: %3.4f |New Expectation: %3.6f\n", change_accepted, bad_change,acceptance_rate,E_old);
		if(acceptance_rate<0.011) poor_acceptance_count++;
		else poor_acceptance_count = 0;
		if(poor_acceptance_count>LIMIT)
		{
			printf("NO PROGRESS FOR %i TEMP CHANGE ITERATIONS, MOVING TO NEXT START STATE\n", LIMIT);
		       	goto end_loop;
		}

		temperature=temperature*EXP_DECAY;
	}
	end_loop: printf("Post-Optimized Expectation:  %f\n", E_best);
	if(PRINTBEST) printf("End array"), print_best_arrays(k_array, j_array, b_array);
	free(k_array), free(j_array), free(b_array),free(k_best), free(j_best), free(b_best), free(k_temp), free(j_temp), free(b_temp),free(psi), free(bonds);
}




void evolve(int *table, unsigned long long int *b,int num_electrons,int N, int lattice[][NX], double *psi, double *k_array, double *j_array, double *b_array, int *bonds)
/*evolve a starting state, psi, by acting on it with exp(ham_dev*-i*time_step). The resulting state is updated as psi and the process repeats TOTAL_STEPS times (TOTAL_TIME/TOTAL_STEPS) until the final psi state is produced. The function contains two methods for calculating exp(ham_dev*-i+time_step), one is a diagonalization method, the other a Pade approximation*/
{
	int i,j;
	double *ham_dev,*ham_t_i, *ham_diag,*exp_matrix,*D, *Vdag;
	ham_dev = (double *) malloc (2*N*N*sizeof (double));
	exp_matrix = (double *) malloc (2*N*N*sizeof (double));
	ham_t_i = (double *) malloc (2*N*N*sizeof (double));
	ham_diag = (double *) malloc(N*N*sizeof(double));
	Vdag = (double*) malloc(N*N*sizeof(double));
	D = (double*) malloc(N*sizeof(double));

	for (i=0; i<TOTAL_STEPS;i++)
	{
		construct_device_hamiltonian(table, b, num_electrons, N, ham_dev, lattice, k_array, j_array, b_array, bonds,i,DEVICE_DIMENSION);
		if(DIAG)
		{
			for (j=0; j<N*N; j++) ham_diag[j]=0.0,Vdag[j]=0.0;
			for (j=0; j<N; j++) D[j]=0.0;
			for (j=0; j<N*N; j++) ham_diag[j] = ham_dev[2*j];//converting an all-real-valued complex matrix into just real matrix
			diag_hermitian_real_double(N, ham_diag,Vdag, D);
			exp_diaganolized_mat(ham_dev, Vdag, D, N);//This function exponentiates D to e^(-iTIME_STEPD)
			matrix_vector_mult(ham_dev,psi, N);
		}
		else
		{
			for (j=0; j<N*N*2; j++) exp_matrix[j] = 0.0, ham_t_i[j]=0.0;
			for (j=0; j<N*N*2; j++) ham_t_i[(j+1)%(N*N*2)] = (ham_dev[j]*-TIME_STEP); //multiplying by -i*dt for the Pade approximation
			exp_general_complex_double(N, ham_t_i, exp_matrix);
			matrix_vector_mult(exp_matrix, psi, N);
		}
	}
	free(ham_dev), free(ham_t_i), free(exp_matrix), free(D), free(Vdag),free(ham_diag);
}




double cost(double *psi, double *ham_mod, int N)
/*Computing the expectation value between ham_mod and psi, <psi|ham_mod|psi>*/
{
	int i=0,INCX = 1,INCY = 1,j;
	double *psi_conj, result=0, resulti=0;
	psi_conj = (double*) malloc (N*2*sizeof(double));
	memcpy(psi_conj, psi, 2*N*sizeof(double));

/*
	double *test1, *test2, sum=0, sumi=0;
	test1 = (double*) malloc(2*N*sizeof(double));
	test2 = (double*) malloc(2*N*sizeof(double));


	for(i=0;i<N;i++) test2[2*i] = psi[2*i];
	for(i=0;i<N;i++) test2[2*i+1] = -psi[2*i+1];//the congugate
	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++){
			sum += ham_mod[2*(N*j+i)]*psi[2*j]-ham_mod[2*(N*j+i)+1]*psi[2*j+1];
			sumi += ham_mod[2*(N*j+i)]*psi[2*j+1]+ham_mod[2*(N*j+i)+1]*psi[2*j];
		}
		test1[2*i]=sum, test1[2*i+1]=sumi;
		sum =0, sumi=0;

	}
	for(i=0;i<N;i++) {
		result+= test1[2*i]*test2[2*i]-test1[2*i+1]*test2[2*i+1];
		resulti+= test1[2*i+1]*test2[2*i]+test1[2*i]*test2[2*i+1];
	}
	if(resulti > 0.00001) printf("  RESOULTTTII %f\n\n\n",resulti);
*/
	result=0;
	matrix_vector_mult(ham_mod, psi, N);//H*psi, the operator acting on the ket, storing in psi
	result = zdotc_(&N, psi_conj, &INCX, psi, &INCY);//psi* * psi, the dot between the complex conj and the result of H*psi
	free(psi_conj);
	return result;
}




double calc_initial_temp(int *table, unsigned long long int *b, int num_electrons, int N, int lattice[][NX], double*ham_mod, gsl_rng *r)
/*Finding an good initial temp that allows an average increase acceptance probability of about ACCEPTANCE_PROB (typically 0.8). Choosing an initial temperature that allows the chance of accepting j,b, and k_array values which increase the expectation value to be ~80%, https://www.phy.ornl.gov/csep/mo/node32.html*/
{
	printf("\n...Calculating initial temperature based on %i random starting states...\n", RANDOM_STATES);
	int *bonds, i=0,j=0, random_row=0,random_col=0, start_state=0, count=0;
	double *psi,*psi_start, *k_array, *j_array,*b_array, E_old=0, E_new=0,sum = 0, change_pm=0, initial_temp=0;
	bonds = (int*) malloc(NUM_SITES*NUM_SITES*sizeof(int));
	k_array = (double *) malloc(2*NUM_SITES*TOTAL_STEPS*sizeof(double));
	j_array = (double *) malloc(2*NUM_SITES*TOTAL_STEPS*sizeof(double));
	b_array = (double *) malloc(NUM_SITES*TOTAL_STEPS*sizeof(double));
	psi_start = (double *) malloc (2*N*sizeof(double));
	psi = (double *) malloc (2*N*sizeof(double));
	assign_bonds(bonds, lattice);
	for (i=0; i<N*2;i++) psi_start[i] =0.0,psi[i]=0.0;
	for (j=0;j<RANDOM_STATES;j++)
	{
		for (i=0; i<N*2;i++) psi_start[i] =0.0;
		start_state = floor(get_random(0,N,r));
		psi_start[start_state] = 1;
		memcpy(psi,psi_start, 2*N*sizeof(double));
		init_arrays(k_array, j_array, b_array,r);

		evolve(table,b,num_electrons,N,lattice, psi, k_array,j_array,b_array, bonds);
		E_old = cost(psi, ham_mod, N);

		for (i=0; i<SWEEPS;i++)
		{
			change_pm = pow(-1,(int)floor(get_random(0,10,r))) * CHANGE;
			random_row = floor(get_random(0,TOTAL_STEPS,r));
			random_col = floor(get_random(0,NUM_SITES*2,r));
			change_array(k_array,j_array,b_array,random_row,random_col,change_pm,i);
			memcpy(psi,psi_start, 2*N*sizeof(double));//resetting psi
			evolve(table,b,num_electrons,N,lattice, psi, k_array,j_array,b_array, bonds);
			E_new = cost(psi, ham_mod, N);

			if (E_new>E_old) sum += (E_new-E_old), count++;
			E_old=E_new;
		}
	}
	free(k_array), free(j_array), free(b_array), free(psi), free(psi_start), free(bonds);
	initial_temp = -(sum/(count*log(ACCEPTANCE_PROB)));
	printf("#######################################################################");
	printf("\nELECTRONS: %i\nDIMENSION: %i\nINIT_TEMP: %f\n",num_electrons, N,initial_temp);
	printf("#######################################################################\n");
	return initial_temp;
}




void construct_device_hamiltonian(int *table, unsigned long long int *b,int num_electrons,int N, double *ham_dev, int lattice[][NX],double *k_array, double *j_array, double *b_array, int *bonds,  int index, int D)
/*constructiing the hamiltonina matrix for the device hamiltonian, using the index-ith row of each j, b, and k array*/
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




void exp_diaganolized_mat(double *ham, double *Vdag, double* D, int N)
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
	for (i =0; i<N; i++) exp_D[2*(i*N+i)] = cos(-TIME_STEP*D[i]), exp_D[2*(i*N+i)+1] = sin(-TIME_STEP*D[i]);
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
	for (i=0;i<N*2;i++) result[i]=0.0;

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
{
	int i, min_i;
	double *D, *Vdag, min_E=1000;
	Vdag = (double*) malloc(N*N*sizeof(double));
	D = (double*) malloc(N*sizeof(double));

	diag_hermitian_real_double(N, ham,Vdag, D);
	for(i=0; i<2*N; i++) ground_state[i] = 0.0;
	for(i=0; i<N; i++) if(D[i]<min_E) min_E = D[i], min_i = i;
	for(i=0; i<N; i++) ground_state[i*2] = Vdag[min_i*N+i];
	free(Vdag), free(D);
}




double get_ground_E(int N, double *ham)
{
	int i;
	double *D, *Vdag, min_E=1000;
	Vdag = (double*) malloc(N*N*sizeof(double));
	D = (double*) malloc(N*sizeof(double));

	diag_hermitian_real_double(N, ham,Vdag, D);
	for(i=0; i<N; i++) if(D[i]<min_E) min_E = D[i];
	free(Vdag), free(D);
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




void copy_arrays(int N, double *k_array, double *j_array, double* b_array,  double* k_to,  double* j_to, double* b_to)
/*storing the updated k, j, and b values*/
{
	memcpy(k_to, k_array, 2*NUM_SITES*TOTAL_STEPS*sizeof(double));
	memcpy(j_to, j_array, 2*NUM_SITES*TOTAL_STEPS*sizeof(double));
	memcpy(b_to, b_array, NUM_SITES*TOTAL_STEPS*sizeof(double));
}




void init_arrays(double *k_array,double *j_array,double *b_array,gsl_rng * r)
/*Initializng the values of the k, j, and b lists which hold the values of the constants for each site (b_array) and between each bond (k_array and j_array)*/
{
	int i,j;
	double upj=1.0,lowj=-1.0,upk=1.0,lowk=-1.0,upb=1.0,lowb=-1.0;
	if(DEVICE_SAME_MODEL) upj=J2,lowj=J1,upk=K2,lowk=K1,upb=B2,lowb=B1;
	for (i=0; i<TOTAL_STEPS;i++)
	{
		for (j=0; j<NUM_SITES*2; j++)
		{
			j_array[i*NUM_SITES*2+j] = fmod((0.25*j+i/0.69),min(abs(lowj),abs(upj)))*pow(-1,j);
			k_array[i*NUM_SITES*2+j] = fmod((0.53*j+i/0.39),min(abs(lowk),abs(upk)))*pow(-1,j+1);
		}
		for (j=0; j<NUM_SITES; j++)
		{
			b_array[i*NUM_SITES+j]= fmod((0.29*j+i/0.69),min(abs(lowb),abs(upb)))*pow(-1,j);
		}
	}
}




void change_array(double *k_array, double *j_array, double *b_array, int random_row, int random_col, double change_pm, int i)
/*changing the j, k, and b arrays. Use the defines at start of program to determine which manipulation functions will be used, and in what order*/
{
	int mod = VARIATIONS;
	if(i%mod==ROW-1) change_row(k_array,j_array,b_array,random_row,change_pm, 1, 0);
	else if(i%mod==COL-1) change_col(k_array,j_array,b_array,random_col,change_pm, 1, 0);
	else if(i%mod==ALTROW-1) change_row(k_array,j_array,b_array,random_row,change_pm, 2, 0);
	else if(i%mod==ALTCOL-1) change_col(k_array,j_array,b_array,random_col,change_pm, 2, 0);
	else if(i%mod==ALTROW2-1) change_row(k_array,j_array,b_array,random_row,change_pm, 2, 1);
	else if(i%mod==ALTCOL2-1) change_col(k_array,j_array,b_array,random_col,change_pm, 2, 1);
	else if(i%mod==SINGLE-1) change_single(k_array,j_array,b_array,random_row,random_col,change_pm);
	else printf("NOPE\n");
}




void change_row(double *k_array,double *j_array,double *b_array, int row, double change, int jump, int offset)
/*Changing all of the lists by a value change at row row. Used in the optimize function*/
{
	int i;
	double upj=1.0,lowj=-1.0,upk=1.0,lowk=-1.0,upb=1.0,lowb=-1.0;
	if(DEVICE_SAME_MODEL) upj=J2,lowj=J1,upk=K2,lowk=K1,upb=B2,lowb=B1;

	for (i=offset; i<NUM_SITES*2; i+=jump)
	{
		if ((lowk < k_array[NUM_SITES*2*row+i] + change) &&  (k_array[NUM_SITES*2*row+i] + change < upk)) k_array[NUM_SITES*2*row+i] += change;
		if ((lowj < j_array[NUM_SITES*2*row+i] + change) && (j_array[NUM_SITES*2*row+i] + change < upj)) j_array[NUM_SITES*2*row+i] += change;
	}
	for (i=offset; i<NUM_SITES; i+=jump)
	{
		if ((lowb < b_array[NUM_SITES*row+i] + change) && (b_array[NUM_SITES*row+i] + change < upb)) b_array[NUM_SITES*row+i] += change;
	}
}




void change_col(double *k_array,double *j_array,double *b_array, int col, double change, int jump, int offset)
/*Changing all of the lists by a value change at col col. Used in the optimize function*/
{
	int i;
	double upj=1.0,lowj=-1.0,upk=1.0,lowk=-1.0,upb=1.0,lowb=-1.0;
	if(DEVICE_SAME_MODEL) upj=J2,lowj=J1,upk=K2,lowk=K1,upb=B2,lowb=B1;

	for (i=offset; i<TOTAL_STEPS; i+=jump)
	{
		if ((lowj < j_array[NUM_SITES*2*i+col] + change) && (j_array[NUM_SITES*2*i+col] + change < upj)) j_array[NUM_SITES*2*i+col] += change;
		if ((lowk < k_array[NUM_SITES*2*i+col] + change) && (k_array[NUM_SITES*2*i+col] + change  < upk)) k_array[NUM_SITES*2*i+col] += change;
	}
	for (i=offset; i<TOTAL_STEPS; i+=jump)
	{
		if ((lowb < b_array[NUM_SITES*i+(int)floor(col/2.0)] + change) && (b_array[NUM_SITES*i+(int)floor(col/2.0)] + change < upb)) b_array[NUM_SITES*i+(int)floor(col/2.0)] += change;
	}
}




void change_single(double *k_array,double *j_array,double *b_array, int row,int col, double change)
/*Changing a single element in each array, at column col, row row*/
{
	double upj=1.0,lowj=-1.0,upk=1.0,lowk=-1.0,upb=1.0,lowb=-1.0;
	if(DEVICE_SAME_MODEL) upj=J2,lowj=J1,upk=K2,lowk=K1,upb=B2,lowb=B1;


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




//double get_random(double lower, double upper)
/*randomly generating a number in [lower, upper), including lower, up to upper (but not including)*/
/*{
	const gsl_rng_type * TT;
	gsl_rng * r;
	struct timeval tv;
	double u, seed;

	gsl_rng_env_setup();
	gettimeofday(&tv,0);
	seed = tv.tv_usec;
	TT = gsl_rng_default;
	r = gsl_rng_alloc (TT);
	gsl_rng_set(r, seed);

	u = gsl_rng_uniform(r);
	gsl_rng_free (r);

	return u*(upper-lower)+lower;
}
*/
double get_random(double lower, double upper, gsl_rng * r)
/*randomly generating a number in [lower, upper), including lower, up to upper (but not including)*/
{
	double u;

	u = gsl_rng_uniform(r);
	return u*(upper-lower)+lower;
}





void print_best_arrays(double *k_array, double *j_array, double *b_array)
{
	int i,j;
	printf("\nK_array:");
	for(i=0;i<TOTAL_STEPS;i++)
	{
		printf("\n");
		for(j=0;j<2*NUM_SITES;j++)
		{
			printf("%5.2f|",k_array[i*NUM_SITES+j]);
//			if(k_array[i*NUM_SITES+j] > 1 || -1 > k_array[i*NUM_SITES+j]) printf("\n\n\n AQJNHDASLJHDAS \n\n\n");
		}
	}
/*	printf("\nJ_array:");
	for(i=0;i<TOTAL_STEPS;i++)
	{
		printf("\n");
		for(j=0;j<2*NUM_SITES;j++)
		{
			printf("%5.2f|",j_array[i*NUM_SITES+j]);
//			if(k_array[i*NUM_SITES+j] > 1 || -1 > k_array[i*NUM_SITES+j]) printf("\n\n\n AQJNHDASLJHDAS \n\n\n");
		}
	}
	printf("\nB_array:");
	for(i=0;i<TOTAL_STEPS;i++)
	{
		printf("\n");
		for(j=0;j<NUM_SITES;j++)
		{
			printf("%5.2f|",b_array[i*NUM_SITES+j]);
//			if(b_array[i*NUM_SITES+j] > 1 || -1 > b_array[i*NUM_SITES+j]) printf("\n\n\n AQJNHDASLJHDAS \n\n\n");
		}
	}*/
	printf("\n");
}




void print_hamiltonianR(double* hamiltonian, int N)
/*prints a real-valued hamiltonian in matrix form*/
{
	int i,j;
	for (i=0;i<N;i++)
	{
		for (j=0;j<N;j++) printf("%4.1f|",(hamiltonian[j*N+i])+0.0);
		printf("\n");
	}
	printf("\n");
}




void print_hamiltonian(double* hamiltonian, int N, bool print_imaginary)
/*prints a complex hamiltonian in matrix form, with the option to just print out the real values*/
{
	int i,j,scalar;
	scalar=2;
	if(print_imaginary) scalar =1;
	for (i=0;i<N*2;i+=scalar)
	{
		if (i%2 == 0)printf("Real: ");
		else printf("Imag: ");
		for (j=0;j<N*2;j+=2) printf("%10.7f|",(hamiltonian[j*N+i]+0.0));
		printf("\n");
	}
	printf("\n");
}




void test_function()
{
	double *matrix, *matrix2, *psi;
	int i, N = 3;
	matrix = (double*) malloc(2*N*N*sizeof(double));
	matrix2 = (double*) malloc(2*N*N*sizeof(double));
	psi = (double*) malloc(2*N*sizeof(double));
	for(i=0;i<N*N*2;i++) matrix[i] = 0.0, matrix2[i] =0.0;
	for(i=0;i<N*2;i++) psi[i] = 0.0;
	matrix[0] = 1;
	matrix[1] = 1;
	matrix[6] = 1;
	matrix[12] = 1;
	matrix[14] = 1;
	matrix[16] = 1;
	psi[1] = 1;
	psi[0] = 1;
	psi[2] = 1;
	psi[4] = 1;
	exp_general_complex_double(N,matrix,matrix2);
	for(i=0;i<2*N;i++) printf("%f\n", psi[i]);
	print_hamiltonian(matrix2, N, true);
	printf("%f\n", cost(psi, matrix, N));
	free(matrix), free(psi);
}
