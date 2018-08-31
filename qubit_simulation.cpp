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

#define OPEN true
#define NX 3 
#define NY NX
#define NUM_SITES NY*NX
#define ELECTRONS 2
#define V 0.2
#define T 0.2
#define TIME_STEP 0.5
#define TOTAL_TIME 10
#define TOTAL_STEPS TOTAL_TIME/TIME_STEP
#define START_STATE 1
#define OPT_ITERATIONS 3
#define CHANGE -0.03
#define TEMP 1
#define PROB_LIMIT 0.5

/*Questions:
negative expectation value?
V-term calculation evaluates open-boundary neighbors?
V and T values?


TODO:
- If >0, keep with prob e^(-deltaE/T)

-carefully check calculations
-double check + vs - in calc
*/


using namespace std;


void optimize(int *table, unsigned long long int *b, int num_electrons, int dimension, int lattice[][NX],double *ham_dev, double*ham_mod);
void evolve(int *table, unsigned long long int *b,int num_electrons,int dimension,  int lattice[][NX], double *ham_dev, double *psi, double *k_list, double *j_list, double *b_list, int *bonds);
double cost(double *psi, double *ham_mod, int dimension);
void construct_device_hamiltonian(int *table, unsigned long long int *b,int num_electrons,int dimension, double *ham_dev, int lattice[][NX],double *k_list, double *j_list, double *b_list, int *bonds, int index);
void construct_model_hamiltonian(int *table, unsigned long long int *b,int num_electrons,int dimension, double *ham_mod,  int lattice[][NX]);
void diag_hermitian_real_double(int N,  double *A, double *Vdag,double *D);
void exp_diaganolized_mat(double *ham_real, double *Vdag, double* D, int dimension);
void exp_general_complex_double(int N, double *A, double *B);
void matrix_vector_mult(double *exp_matrix, double *psi, int dimension);
int hop(unsigned long long int b, unsigned long long int *v,int n, int j);
int find(int dimension,unsigned long long int *v, unsigned long long int *b);
unsigned long long int choose(int num_electrons);
int combinations ( int num_electrons, unsigned long long int *b,int *tab, int dimension);
void assign_bonds(int *bonds, int lattice[][NX]);
void init_lists(double *k_list, double *j_list,double *b_list);
void change_lists(double *k_list, double *j_list,double *b_list, int index, double change);
void construct_lattice(int lattice[][NX]);
int get_neighbors(int site, int *neighbors, int lattice[][NX]);
int get_random(int lower, int upper);
void print_hamiltonian(double* hamiltonian, int dimension, char print_imaginary);
void print_hamiltonianR(double *hamiltonian, int dimension);
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
	int *table,num_electrons,dimension,lattice[NY][NX];
	unsigned long long int *b;
	double *ham_mod,*ham_dev;
	num_electrons=ELECTRONS;
	dimension=choose(num_electrons);

	b = (unsigned long long int*) malloc(dimension*sizeof(double));
	table=(int*) malloc(num_electrons*dimension*sizeof(int));
	ham_mod = (double *) malloc (2*dimension*dimension*sizeof (double));
	ham_dev = (double *) malloc (2*dimension*dimension*sizeof (double));

	construct_lattice(lattice);
	combinations (num_electrons,b,table, dimension);
	construct_model_hamiltonian(table, b, num_electrons, dimension, ham_mod, lattice);

	optimize(table, b, num_electrons, dimension, lattice, ham_dev, ham_mod);

	free(b), free(table), free(ham_mod), free(ham_dev);
	exit (0);
}





void optimize(int *table, unsigned long long int *b, int num_electrons, int dimension, int lattice[][NX],double *ham_dev, double*ham_mod)
/*Optimize the values in the j, b, and k list in order to produce the lowest energy (expectation value) between the final state (psi), produced by evolve, and the model hamiltonian. This is done by randomly selecting one row of each list, making a slight change, then determining if the new energy is lower than the old. If the new energy is greater than the old, keep with probability exp(delta_E/Temp)*/
{
	int i,j, random_int, *bonds;
	double *psi,*psi_reset, E_old, E_new, delta_E,*k_list, *j_list,*b_list;
	bonds = (int*) malloc(NUM_SITES*NUM_SITES*sizeof(int));
	k_list = (double *) malloc(2*NUM_SITES*TOTAL_STEPS*sizeof(double));
	j_list = (double *) malloc(2*NUM_SITES*TOTAL_STEPS*sizeof(double));
	b_list = (double *) malloc(NUM_SITES*TOTAL_STEPS*sizeof(double));
	psi_reset = (double *) malloc (2*dimension*sizeof(double));
	psi = (double *) malloc (2*dimension*sizeof(double));

	for (i=0; i<dimension*2;i++) psi_reset[i] =0.0;
	psi_reset[START_STATE] = 1;
	memcpy(psi,psi_reset, 2*dimension*sizeof(double));

	assign_bonds(bonds, lattice);
	init_lists(k_list, j_list, b_list);

	evolve(table,b,num_electrons,dimension,lattice, ham_dev, psi, k_list,j_list,b_list, bonds);
	E_old = cost(psi, ham_mod, dimension);
	printf("\nExpectation: %f\n\n", E_old);

	for (i=0; i<OPT_ITERATIONS;i++)
	{
		random_int = get_random(0,TOTAL_STEPS-1);	
		change_lists(k_list,j_list,b_list,random_int,CHANGE);
		
		memcpy(psi,psi_reset, 2*dimension*sizeof(double));//resetting psi
		evolve(table,b,num_electrons,dimension,lattice, ham_dev, psi, k_list,j_list,b_list, bonds);
		E_new = cost(psi, ham_mod, dimension);
		delta_E = E_new-E_old;

		if (E_new<E_old) E_old=E_new;
		else if (exp(-delta_E/TEMP)>PROB_LIMIT) change_lists(k_list,j_list,b_list,random_int,-CHANGE), printf("exp(x)>Lim\n");
			
	}	
	printf("\nExpectation: %f\n\n", E_old);
	free(k_list), free(j_list), free(b_list), free(psi), free(psi_reset), free(bonds);
}




void evolve(int *table, unsigned long long int *b,int num_electrons,int dimension, int lattice[][NX], double *ham_dev, double *psi, double *k_list, double *j_list, double *b_list, int *bonds)
/*evolve a starting state, psi, by acting on it with exp(ham_dev*-i*time_step). The resulting state is updated as psi and the process repeats TOTAL_STEPS times (TOTAL_TIME/TOTAL_STEPS) until the final psi state is produced. The function contains two methods for calculating exp(ham_dev*-i+time_step), one is a diagonalization method, the other a Pade approximation*/
{
	int i,j;
	double *ham_t_i, *ham_diag,*exp_matrix,*D, *Vdag;
	ham_diag = (double *) malloc(dimension*dimension*sizeof(double));
	D = (double*) malloc(sizeof(double)*dimension);
	Vdag = (double*) malloc(sizeof(double)*dimension*dimension);
	ham_t_i = (double *) malloc (2*dimension*dimension*sizeof (double));
	exp_matrix = (double *) malloc (2*dimension*dimension*sizeof (double));


	for (i=0; i<TOTAL_STEPS;i++)
	{

		construct_device_hamiltonian(table, b, num_electrons, dimension, ham_dev, lattice, k_list, j_list, b_list, bonds,i); 

		for (j=0; j<dimension*dimension*2; j++) exp_matrix[j] = 0.0,exp_matrix[j]=0.0;
		for (j=0; j<dimension*dimension; j++) ham_diag[j]=0.0,Vdag[j]=0.0;
		for (j=0; j<dimension; j++) D[j]=0.0;
		for (j=0; j<dimension*dimension*2; j++) ham_t_i[(j+1)%(dimension*dimension*2)] = (ham_dev[j]*-TIME_STEP)+0.0; //multiplying by -i*dt for the Pade approximation
		for (j =0; j<dimension*dimension; j++) ham_diag[j] = ham_dev[2*j];//converting an all-real-valued complex matrix into just real matrix

		//The diagonalization method
		//diag_hermitian_real_double(dimension, ham_diag,Vdag, D);
		//exp_diaganolized_mat(ham_dev, Vdag, D, dimension);//This function exponentiates D to e^(-iTIME_STEPD)
		//matrix_vector_mult(ham_dev,psi, dimension);


		//The Pade approximation
		exp_general_complex_double(dimension, ham_t_i, exp_matrix);
		matrix_vector_mult(exp_matrix, psi, dimension);

	}
	free(ham_diag), free(D), free(Vdag), free(ham_t_i), free(exp_matrix);
}




double cost(double *psi, double *ham_mod, int dimension)
/*Computing the expectation value between ham_mod and psi, <psi|ham_mod|psi>*/
{
	int i;
	double *psi_conj, result;
	int INCX = 1;
	int INCY = 1;
	psi_conj = (double*) malloc (dimension*2*sizeof(double));
	memcpy(psi_conj, psi, 2*dimension*sizeof(double));

	matrix_vector_mult(ham_mod, psi, dimension);//H*psi, the operator acting on the ket
	result = zdotc_(&dimension, psi_conj, &INCX, psi, &INCY);//psi* * psi, the dot between the complex conj and the result of Hpsi
	free(psi_conj);
	return result;
}




void construct_device_hamiltonian(int *table, unsigned long long int *b,int num_electrons,int dimension, double *ham_dev, int lattice[][NX],double *k_list, double *j_list, double *b_list, int *bonds,  int index)
/*constructiing the hamiltonina matrix for the device hamiltonian, using the index-ith row of each j, b, and k array*/
{
	int i,ii,j,x,y,state,site,sign,bond,neighbor_count,*neighbors;
	unsigned long long int *v,comparison;
	v = (unsigned long long int*) malloc(sizeof(unsigned long long int));
	neighbors = (int*) malloc(4*sizeof(int));
	for (i=0; i<dimension*dimension*2; i++) ham_dev[i] = 0;

	for (i=0;i<dimension;i++)
	{
		for (j=0;j<num_electrons;j++)//The J term calculation
		{
			site=table[i*num_electrons+j];
			neighbor_count = get_neighbors(site, neighbors, lattice);
			for (ii=0; ii<neighbor_count; ii++)
			{
				if (((1ULL<<(neighbors[ii]-1))&b[i])==0)//making sure neighbor is not occupied, otherwise nothing happens
				{
					hop(b[i], v,site, neighbors[ii]);
					state=find(dimension,v, b);
					bond = bonds[(site-1)*NUM_SITES+neighbors[ii]-1];
					ham_dev[(dimension*i+state-1)*2] -= j_list[index*NUM_SITES*2+bond];
				}
			}
		}
		for (j=1;j<(NUM_SITES);j++)//The K term calculation
		{
			site=j;
			neighbor_count = get_neighbors(site, neighbors, lattice);
			for (ii=0; ii<neighbor_count;ii++)
			{
				if (neighbors[ii] > site)
				{
					sign = -1;
					comparison = (1ULL<<(neighbors[ii]-1))+(1ULL<<(site-1));
					if((comparison&b[i])==comparison || (comparison&b[i])==0) sign = 1;
					bond = bonds[(site-1)*NUM_SITES+neighbors[ii]-1];
					ham_dev[((dimension*i)+i)*2] += k_list[index*NUM_SITES*2+bond]*sign;
				}
			}
		}
		for (j=0; j<NUM_SITES;j++)//The B term calculation, DOUBLE CHECK THIS, changing eval of if changes result
		{
			sign = -1;
			if(((1ULL<<j)&b[i])>0) sign=1;
			ham_dev[((dimension*i)+i)*2] += b_list[index*NUM_SITES+j]*sign;
		}
	}
	free(neighbors);
}





void construct_model_hamiltonian(int *table, unsigned long long int *b,int num_electrons,int dimension, double *ham_mod, int lattice[][NX])
/*Constructing the hamiltonian matrix for the model hamiltonian*/
{
	int i,ii,j,x,y,state,site,neighbor,sign,neighbor_count, *neighbors;
	unsigned long long int *v,comparison;
	v = (unsigned long long int*) malloc(sizeof(unsigned long long int));
	neighbors = (int*) malloc(sizeof(int)*4);
	for (i=0; i<dimension*dimension*2; i++) ham_mod[i] = 0;

	for (i=0;i<dimension;i++)
	{
		for (j=0;j<num_electrons;j++)//The T term calculation
		{
			site=table[i*num_electrons+j];

			neighbor_count = get_neighbors(site, neighbors, lattice);
			for (ii=0; ii<neighbor_count; ii++)
			{
				if (((1ULL<<(neighbors[ii]-1))&b[i])==0)//checking if the neighbor is occupied
				{
					sign = hop(b[i], v,site, neighbors[ii]);
					if (sign==0) sign=1;
					else sign=-1;
					state=find(dimension,v, b);
					ham_mod[(dimension*i+state-1)*2] -= (T*sign);
				}
			}
		}

		for (j=1;j<(NUM_SITES);j++)//The V term calculation

		{
			site=j;
			neighbor_count = get_neighbors(site, neighbors, lattice);
			for (ii=0; ii<neighbor_count;ii++)
			{
				if (neighbors[ii] > site)
				{
					sign = -1;
					comparison = (1ULL<<(neighbors[ii]-1))+(1ULL<<(site-1));
					if((comparison&b[i])==comparison || (comparison&b[i])==0) sign = 1;
					ham_mod[((dimension*i)+i)*2] += sign*V;
				}
			}
		}
	}
	free(neighbors);
}





void diag_hermitian_real_double(int N,  double *A, double *Vdag,double *D)
/*diagonalizing an real square matrix A. Stores the eigenvectors in Vdag and the eigenvalues in D*/
{
	char JOBZ,UPLO;
	int LDA,LWORK, INFO;
	double *WORK;
	JOBZ='V';
	UPLO='U'; 
	LDA=N; 
	LWORK=-1; 
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




void exp_diaganolized_mat(double *ham, double *Vdag, double* D, int dimension)
/*calculating the exponential of a diagonalized decomposition where the matrix A = Vdag*D*Vdag_inv. calculating Vdag*exp(D)*Vdag_inv =exp(A), storing the result in ham*/
{
	int i, *IPIV, N,LWORK, INFO, LDA, LDB, LDC;
	double *exp_D, *temp_mat, *Vdag_z, *Vdag_z_inv, *WORK;
	N = dimension;
	LWORK = N*N;
	LDA=N;
	LDB=N;
	LDC=N;
	char TRANSA = 'N';
	char TRANSB = 'N';
	double ALPHA[2];
	ALPHA[0]=1.0;
	ALPHA[1]=0.0;
	double BETA[2];
	BETA[0]=0.0;
	BETA[1]=0.0;

	IPIV = (int*) malloc(N*sizeof(int));
	WORK = (double*) malloc(LWORK*sizeof(double));
	exp_D = (double *) malloc (2*dimension*dimension*sizeof (double));
	temp_mat = (double *) malloc (2*dimension*dimension*sizeof (double));
	Vdag_z = (double *) malloc (2*dimension*dimension*sizeof (double));
	Vdag_z_inv = (double *) malloc (2*dimension*dimension*sizeof (double));

	for (i =0; i<dimension*dimension*2; i++) exp_D[i] = 0, temp_mat[i] = 0,Vdag_z[i]=0, Vdag_z_inv[i]=0, ham[i]=0;
	for (i =0; i<dimension*dimension; i++) Vdag_z[2*i] = Vdag[i];
	for (i =0; i<dimension; i++) exp_D[2*(i*dimension+i)] = cos(-TIME_STEP*D[i]), exp_D[2*(i*dimension+i)+1] = sin(-TIME_STEP*D[i]);
	dgetrf_(&N,&N,Vdag,&N,IPIV,&INFO);
	dgetri_(&N,Vdag,&N,IPIV,WORK,&LWORK,&INFO);//inverting vdag
	for (i =0; i<dimension*dimension; i++) Vdag_z_inv[2*i] = Vdag[i];

	zgemm_(&TRANSA, &TRANSB, &N, &N, &N, ALPHA, Vdag_z, &LDA, exp_D, &LDB, BETA, temp_mat, &LDC); //matrix mult
	zgemm_(&TRANSA, &TRANSB, &N, &N, &N, ALPHA, temp_mat, &LDA, Vdag_z_inv, &LDB, BETA, ham, &LDC); //matrix mult
	free(exp_D),free(temp_mat),free(Vdag_z), free(Vdag_z_inv), free(IPIV), free(WORK);
}




void exp_general_complex_double(int N, double *A, double *B)
/*Using the Pade Approximation to calulate exp(A), where A is a complex matrix storing real values at even indices and their imaginary parts at index +1. Storing the result in B*/
{
	int M,K,ii,jj,kk,s,p,q,INFO,LDA,LDB,LDC,NRHS, *IPIV;
	double *row_norm,*X,*Y,*Z,*E,*D;
	double norm,c;
	char TRANSA='N';
	char TRANSB='N';
	char TRANS='N';

	M=N;
	K=N;
	double ALPHA[2];
	ALPHA[0]=1.0;
	ALPHA[1]=0.0;
	double BETA[2];
	BETA[0]=0.0;
	BETA[1]=0.0;
	LDA=N;
	LDB=N;
	LDC=N;
	NRHS=N;

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
	c=0.5;
	p=1;
	q=6;


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



void matrix_vector_mult(double *matrix, double *psi, int dimension)
/*Matrix vector multiplication, psi=matrix*psi*/
{
	int i;
	double *result;
	char TRANS = 'N';
	int M = dimension;
	int N = dimension;
	int LDA = dimension;
	int INCX = 1;
	int INCY = 1;
	double ALPHA[2];
	ALPHA[0]=1.0;
	ALPHA[1]=0.0;
	double BETA[2];
	BETA[0]=0.0;
	BETA[1]=0.0;
	result = (double *) malloc (dimension*2*sizeof(double));
	for (i=0;i<dimension*2;i++) result[i]=0.0;

	zgemv_(&TRANS, &M, &N,ALPHA,matrix, &LDA, psi, &INCX, BETA, result, &INCY);
	memcpy(psi,result, 2*dimension*sizeof(double));
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




int find(int dimension,unsigned long long int *v,unsigned long long int *b)
/*find the position of a given combination v (vector of length k) in the table of all combinations tab*/
{
	int first, last, mid;
	first=0;
	last=dimension-1;

	while (first <= last)
	{
		mid = (int) ((first + last) / 2.0);
		if (*v > b[mid]) first = mid + 1;
		else if (*v < b[mid]) last = mid - 1;
		else return mid+1;
	}
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




int combinations ( int num_electrons,  unsigned long long int *b,int *tab, int dimension)
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
	for (i=1;i<dimension;i++)
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
		neighbor_count = get_neighbors(site, neighbors, lattice);
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




void init_lists(double *k_list,double *j_list,double *b_list)
/*Initializng the values of the k, j, and b lists which hold the values of the constants for each site (b_list) and between each bond (k_list and j_list)*/
{
	int i,j;
	for (i=0; i<TOTAL_STEPS;i++)
	{
		for (j=0; j<NUM_SITES*2; j++)
		{
			j_list[i*NUM_SITES*2+j] = 0.05*j*pow(-1,j)+i/20;
			k_list[i*NUM_SITES*2+j] = 0.03*j*pow(-1,j+1)+i/20;;
		}
		for (j=0; j<NUM_SITES; j++)
		{
			b_list[i*NUM_SITES+j]= 0.09*j*pow(-1,j)+i/20;
		}
	}
}




void change_lists(double *k_list,double *j_list,double *b_list, int index, double change)
/*Changing all of the lists by a value change at row index. Used in the optimize function*/
{
	int i;
	for (i=0; i<NUM_SITES*2; i++) j_list[NUM_SITES*2*index+i] += change, k_list[NUM_SITES*2*index+i] +=change;
	for (i=0; i<NUM_SITES; i++) b_list[NUM_SITES*index+i] += change;
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





int get_neighbors(int site, int *neighbors, int lattice[][NX])
/*Gets the neighbors of a site, returning all 4 neighbors if open_boundry coditions are true, otherwise just closed-boudnary neighbors*/
{
	int x,y,i,count=0;
	for (x=0;x<NX;x++) for (y=0;y<NX;y++) if(site == lattice[x][y]) goto end_loop;//Finding the coordinates
	end_loop: for(i=0;i<4;i++) neighbors[i]=0;

	if (OPEN) 
	{
		neighbors[0] = lattice[(x+1)%NX][y]; 
		neighbors[1] = lattice[(x+(NX-1))%NX][y];
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




int get_random(int lower, int upper)
/*randomly generating a number in [lower, upper]*/
{
	const gsl_rng_type * TT;
	gsl_rng * r;
	int i;
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

	return floor(((upper+1)-lower)*u+lower);
}




void print_hamiltonianR(double* hamiltonian, int dimension)
/*prints a real-valued hamiltonian in matrix form*/
{
	int i,j;
	for (i=0;i<dimension;i++)
	{
		for (j=0;j<dimension;j++) printf("%4.1fi",(hamiltonian[j*dimension+i])+0.0);
		printf("\n");
	}
	printf("\n");
}




void print_hamiltonian(double* hamiltonian, int dimension, bool print_imaginary)
/*prints a complex hamiltonian in matrix form, with the option to just print out the real values*/
{
	int i,j,scalar;
	scalar=2;
	if(print_imaginary) scalar =1;
	for (i=0;i<dimension*2;i+=scalar)
	{
		if (i%2 == 0)printf("Real: ");
		else printf("Imag: ");
		for (j=0;j<dimension*2;j+=2) printf("%5.2f|",(hamiltonian[j*dimension+i]+0.0));
		printf("\n");
	}
	printf("\n");
}
