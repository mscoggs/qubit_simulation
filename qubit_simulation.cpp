#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <complex>
#include <math.h>
#include <time.h>

#define NEIGHBORS 4
#define nx 3
#define ELECTRONS 1
#define V .4
#define T .3
#define dT 0.5
#define Tao 10
#define PsiStartState 2
#define NCOUNT Tao/dT


/*Questions:
negative expectation value?


TODO:
GET SSH TO WORK
-change j,k,b list
-get random interval, change slightly, calc expecation
-if delta-E<0, good, keep change. If >0, keep with prob e^(-deltaE/T) 

Less important TODO
-documentation
-algorithmic way to find bonds
-open boundary conditions
-clean up
*/


using namespace std;
typedef struct
{
	int coordX, coordY;
} Coordinates;

int find_bond[][nx*nx]=
	{{0,1,3,0,0,10,16,0,0},
	{1,0,2,0,11,0,0,17,0},
	{3,2,0,12,0,0,0,0,18},
	{0,0,12,0,5,6,0,0,15},
	{0,11,0,5,0,4,0,14,0},
	{10,0,0,6,4,0,13,0,0},
	{16,0,0,0,0,13,0,7,9},
	{0,17,0,0,14,0,7,0,8},
	{0,0,18,15,0,0,9,8,0}};



void evolve(int *table, unsigned long long int *b,int electrons,int dimension, int sites,  int lattice[][nx], double *ham_dev,double *psi);
double cost(double *psi, double *ham_mod, int dimension);
void construct_device_hamiltonian(int *table, unsigned long long int *b,int electrons,int dimension, double *ham_dev, int sites, int lattice[][nx],double *k_list, double *j_list, double *b_list, int index);
void construct_model_hamiltonian(int *table, unsigned long long int *b,int electrons,int dimension, double *ham_mod, int sites,  int lattice[][nx]);
void diag_hermitian_real_double(int N,  double *A, double *Vdag,double *D);
void exp_diaganolized_mat(double *ham_real, double *Vdag, double* D, int dimension);
void exp_general_complex_double(int N, double *A, double *B);
void matrix_vector_mult(double *exp_matrix, double *psi, int dimension);
int hop(unsigned long long int *v1, unsigned long long int *v2,int n, int j);
int find(int dimension,unsigned long long int *v, unsigned long long int *b);
unsigned long long int choose(int n, int k);
int combinations ( int n, int k, unsigned long long int *b,int *tab);
void init_lists(double *k_list, double *j_list,double *b_list);
void change_lists(double *k_list, double *j_list,double *b_list, int index);
void construct_lattice(int lattice[][nx]);
Coordinates find_coordinates(int site_num, int lattice[][nx]);
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
	int i,*table,*v,*v2,sites,electrons,dimension,lattice[nx][nx],*j_list,*b_list,*k_list;
	unsigned long long int *b;
	double *psi,*ham_mod,*ham_dev,expectation;
	sites=nx*nx;
	electrons=ELECTRONS;
	dimension=choose(sites,electrons);

	b = (unsigned long long int*) malloc(dimension*sizeof(double));
	table=(int*) malloc(electrons*dimension*sizeof(int));
	ham_mod = (double *) malloc (2*dimension*dimension*sizeof (double));
	ham_dev = (double *) malloc (2*dimension*dimension*sizeof (double));
	psi = (double *) malloc (2*dimension*sizeof(double));
	for (i=0; i<dimension*2;i++) psi[i] =0.0;
	psi[PsiStartState] = 1;

	construct_lattice(lattice);
	combinations (sites,electrons,b,table);
	construct_model_hamiltonian(table, b, electrons, dimension, ham_mod, sites, lattice);


	evolve(table,b,electrons,dimension,sites,lattice, ham_dev, psi);
	expectation = cost(psi, ham_mod, dimension);
	printf("\nExpectation %f\n\n", expectation);
	free(b), free(table), free(ham_mod), free(ham_dev), free(psi);
	exit (0);
}




void evolve(int *table, unsigned long long int *b,int electrons,int dimension, int sites,  int lattice[][nx], double *ham_dev, double *psi)
{
	int i,j;
	double *ham_t_i, *ham_diag,*exp_matrix,*D, *Vdag,*b_list,*k_list,*j_list;
	k_list = (double *) malloc(2*nx*nx*NCOUNT*sizeof(double));
	j_list = (double *) malloc(2*nx*nx*NCOUNT*sizeof(double));
	b_list = (double *) malloc(nx*nx*NCOUNT*sizeof(double));
	ham_diag = (double *) malloc(dimension*dimension*sizeof(double));
	D = (double*) malloc(sizeof(double)*dimension);
	Vdag = (double*) malloc(sizeof(double)*dimension*dimension);
	ham_t_i = (double *) malloc (2*dimension*dimension*sizeof (double));
	exp_matrix = (double *) malloc (2*dimension*dimension*sizeof (double));

	init_lists(k_list, j_list, b_list);

	for (i=1; i<NCOUNT;i++)
	{

		change_lists(k_list,j_list,b_list,i);
		construct_device_hamiltonian(table, b, electrons, dimension, ham_dev, sites,lattice, k_list, j_list, b_list);

		for (j=0; j<dimension*dimension*2; j++) exp_matrix[j] = 0.0,exp_matrix[j]=0.0;
		for (j=0; j<dimension*dimension; j++) ham_diag[j]=0.0,Vdag[j]=0.0;
		for (j=0; j<dimension; j++) D[j]=0.0;
		for (j=0; j<dimension*dimension*2; j++) ham_t_i[(j+1)%(dimension*dimension*2)] = (ham_dev[j]*-dT)+0.0; //multiplying by -i*dt for the Pade approximation
		for (j =0; j<dimension*dimension; j++) ham_diag[j] = ham_dev[2*j];//converting an all-real-valued complex matrix into just real matrix


		/* diag_hermitian_real_double(dimension, ham_diag,Vdag, D);//The diagonalization method
		exp_diaganolized_mat(ham_dev, Vdag, D, dimension);//This function exponentiates D to e^(-idTD)
		matrix_vector_mult(ham_dev,psi, dimension);*/



		exp_general_complex_double(dimension, ham_t_i, exp_matrix);//The Pade approximation
		matrix_vector_mult(exp_matrix, psi, dimension);

	}
	free(k_list), free(j_list), free(b_list), free(ham_diag), free(D), free(Vdag), free(ham_t_i), free(exp_matrix);
}




double cost(double *psi, double *ham_mod, int dimension)
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




void construct_device_hamiltonian(int *table, unsigned long long int *b,int electrons,int dimension, double *ham_dev, int sites, int lattice[][nx],double *k_list, double *j_list, double *b_list, int index)
/*Constructing the hamiltonian matrix for the device hamiltonian*/
{
	int i,ii,j,x,y,state,site,neighbor,sign,bond;
	unsigned long long int *v1,*v2,comparison;
	Coordinates coord;
	v1=(unsigned long long int*) malloc(sizeof(unsigned long long int));
	v2=(unsigned long long int*) malloc(sizeof(unsigned long long int));

	for (i=0; i<dimension*dimension*2; i++) ham_dev[i] = 0;

	for (i=0;i<dimension;i++)
	{
		for (j=0;j<electrons;j++)//The J term calculation
		{
			site=table[i*electrons+j];
			coord = find_coordinates(site, lattice);
			x = coord.coordX;
			y = coord.coordY;
			int neighbors[NEIGHBORS] = {lattice[(x+1)%nx][y], lattice[(x+(nx-1))%nx][y], lattice[x][(y+1)%nx], lattice[x][(y+(nx-1))%nx]};

			for (ii=0; ii<NEIGHBORS; ii++)
			{
				if (((1ULL<<(neighbors[ii]-1))&b[i])==0)//making sure neighbor is not occupied, otherwise nothing happens
				{
					memcpy(&v1,&b[i], sizeof(unsigned long long int));
					hop(v1, v2,site, neighbors[ii]);
					state=find(dimension,v2, b);
					bond = find_bond[site-1][neighbors[ii]-1];
					ham_dev[(dimension*i+state-1)*2] += j_list[index*nx*nx*2+bond];
				}
			}
		}
		for (j=1;j<(nx*nx);j++)//The K term calculation
		{
			site=j;
			coord = find_coordinates(site, lattice);
			x = coord.coordX;
			y = coord.coordY;
			int neighbors[NEIGHBORS] = {lattice[(x+1)%nx][y], lattice[(x+(nx-1))%nx][y], lattice[x][(y+1)%nx], lattice[x][(y+(nx-1))%nx]};

			for (ii=0; ii<NEIGHBORS;ii++)
			{
				if (neighbors[ii] > site)
				{
					sign = -1;
					comparison = (1ULL<<(neighbors[ii]-1))+(1ULL<<(site-1));
					if((comparison&b[i])==comparison || (comparison&b[i])==0) sign = 1;
					bond = find_bond[site][neighbors[ii]];
					ham_dev[((dimension*i)+i)*2] += k_list[index*nx*nx*2+bond];
				}
			}
		}
		for (j=0; j<nx*nx;j++)//The B term calculation
		{
			if((1ULL<<j)&b[i]) ham_dev[((dimension*i)+i)*2]+=b_list[index*nx*nx+j];
			else ham_dev[((dimension*i)+i)*2]-=b_list[index*nx*nx+j];
		}
	}
}





void construct_model_hamiltonian(int *table, unsigned long long int *b,int electrons,int dimension, double *ham_mod, int sites,  int lattice[][nx])
/*Constructing the hamiltonian matrix for the model hamiltonian*/
{
	int i,ii,j,x,y,state,site,neighbor,sign;
	unsigned long long int *v1,*v2,comparison;
	Coordinates coord;
	v1=(unsigned long long int*) malloc(sizeof(unsigned long long int));
	v2=(unsigned long long int*) malloc(sizeof(unsigned long long int));
	for (i=0; i<dimension*dimension*2; i++) ham_mod[i] = 0;

	for (i=0;i<dimension;i++)
	{
		for (j=0;j<electrons;j++)//The T term calculation
		{
			site=table[i*electrons+j];
			coord = find_coordinates(site, lattice);
			x = coord.coordX;
			y = coord.coordY;
			int neighbors[NEIGHBORS] = {lattice[(x+1)%nx][y], lattice[(x+(nx-1))%nx][y], lattice[x][(y+1)%nx], lattice[x][(y+(nx-1))%nx]};

			for (ii=0; ii<NEIGHBORS; ii++)
			{
				if (((1ULL<<(neighbors[ii]-1))&b[i])==0)//checking if the neighbor is occupied
				{
					memcpy(&v1,&b[i], sizeof(unsigned long long int));
					sign = hop(v1, v2,site, neighbors[ii]);
					if (sign==0) sign=1;
					else sign=-1;
					state=find(dimension,v2, b);
					ham_mod[(dimension*i+state-1)*2] -= (T*sign);
				}
			}
		}
		for (j=1;j<(nx*nx);j++)//The V term calculation
		{
			site=j;
			coord = find_coordinates(site, lattice);
			x = coord.coordX;
			y = coord.coordY;
			int neighbors[NEIGHBORS] = {lattice[(x+1)%nx][y], lattice[(x+(nx-1))%nx][y], lattice[x][(y+1)%nx], lattice[x][(y+(nx-1))%nx]};

			for (ii=0; ii<NEIGHBORS;ii++)
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
}





void diag_hermitian_real_double(int N,  double *A, double *Vdag,double *D)
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
	for (i =0; i<dimension; i++) exp_D[2*(i*dimension+i)] = cos(-dT*D[i]), exp_D[2*(i*dimension+i)+1] = sin(-dT*D[i]);
	dgetrf_(&N,&N,Vdag,&N,IPIV,&INFO);
	dgetri_(&N,Vdag,&N,IPIV,WORK,&LWORK,&INFO);//inverting vdag
	for (i =0; i<dimension*dimension; i++) Vdag_z_inv[2*i] = Vdag[i];

	zgemm_(&TRANSA, &TRANSB, &N, &N, &N, ALPHA, Vdag_z, &LDA, exp_D, &LDB, BETA, temp_mat, &LDC); //matrix mult
	zgemm_(&TRANSA, &TRANSB, &N, &N, &N, ALPHA, temp_mat, &LDA, Vdag_z_inv, &LDB, BETA, ham, &LDC); //matrix mult
	free(exp_D),free(temp_mat),free(Vdag_z), free(Vdag_z_inv), free(IPIV), free(WORK);
}




void exp_general_complex_double(int N, double *A, double *B)
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




int hop(unsigned long long int *v1, unsigned long long int *v2,int n, int j)
/*given a state v1 generates the state v2 obtained by hopping from site v1[n] to site j (which should not be in v1), and outputs the fermionic sign*/
{
	unsigned long long int i,x,y,vtemp;
	int z_count = 0;
	vtemp = (unsigned long long int) v1;
	x = (1ULL << (n-1)) + (1ULL << (j-1));
	for (i=n;i<j-1;i++)  if((1ULL<<i) & (vtemp)) z_count++;
	if((1ULL<<4) & (1ULL<<4));
	y = (x ^ vtemp);
	memcpy(v2, &y, sizeof(unsigned long long int));
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




unsigned long long int choose(int n, int k)
/*Returning the number of combinations of k out of n. (n-choose-k), finding number of unordered combinations with n electrons in k sites.*/
{
	int i;
	unsigned long long int c;
	c=1ULL;
	for (i=0;i<k;i++) c=c*(n-i);
	for (i=0;i<k;i++) c=c/(i+1);
	return c;
}




int combinations ( int n, int k, unsigned long long int *b,int *tab)
/*Returning the number of combinations of k out of n. (n-choose-k), finding number of unordered combinations with n electrons in k sites.*/
{
	unsigned long long int x,y;
	int i,c,d;
	x=0ULL;
	for (i=0;i<k;i++) x=x+(1ULL<<i);
	b[0]=x;
	c=0;
	d=0;
	i=0;
	while ((c<n)&& (d<k))
	{
		if (x & (1ULL<<c))
		{
			tab[i*k+d]=c+1;
			d++;
		}
		c++;
	}
	for (i=1;i<choose(n,k);i++)
	{
		y = (x | (x - 1)) + 1;
		x = y | ((((y & -y) / (x & -x)) >> 1) - 1);
		b[i]=x;
		c=0;
		d=0;
		while ((c<n)&& (d<k))
		{
			if (x & (1ULL<<c))
			{
				tab[i*k+d]=c+1;
				d++;
			}
			c++;
		}
	}
}




void init_lists(double *k_list,double *j_list,double *b_list)
{
	int j;
	for (i=0; i<NCOUNT;i++)
	{		
		for (j=0; j<nx*nx*2; j++)
		{
			j_list[i*nx*nx*2+j] = 0.05*j*pow(-1,j)+i/20;
			k_listt[i*nx*nx*2+j] = 0.03*j*pow(-1,j+1)+i/20;;
		}
		for (j=0; j<nx*nx; j++)
		{
			b_list[i*nx*nx+j]= 0.09*j*pow(-1,j)+i/20;
		}
	}
}	




void change_lists(double *k_list,double *j_list,double *b_list, int index)
{
	if (index==0) return;
	int i;
	for (i=0; i<nx*nx*2; i++) j_list[i] = j_list[i]*.9, k_list[i] = k_list[i]*.5;
	for (i=0; i<nx*nx; i++) b_list[i] = b_list[i]*.9*pow(-1,i);
}




void construct_lattice(int lattice[][nx])
/*Constructing the lattice of dimension nx*nx*/
{
	int x,y;
	for (x=0;x<nx;x++)
	{
		if (x%2 ==1) for (y=0;y<nx;y++) lattice[x][y] = (nx*x)+y+1;
		else for (y=0;y<nx;y++) lattice[x][nx-1-y] = nx*x+y+1;
	}
}




Coordinates find_coordinates(int site_num, int lattice[][nx])
/*Finding the coordinates of a given site*/
{
	int x,y;
	for (x=0;x<nx;x++)
	{
		for (y=0;y<nx;y++)
		{
			if(site_num == lattice[x][y])
			{
				Coordinates coords = { x,y };
				return coords;
			}
		}
	}
}
















//functions for debugging purposes
void print_hamiltonianR(double* hamiltonian, int dimension)
{
	int i,j;
	for (i=0;i<dimension;i++)
	{
		for (j=0;j<dimension;j++) printf("%4.1fi",(hamiltonian[j*dimension+i])+0.0);
		printf("\n");
	}
	printf("\n");
}

void print_hamiltonian(double* hamiltonian, int dimension, char print_imaginary)
{
	int i,j,scalar;
	scalar=2;
	if(print_imaginary =='y') scalar =1;
	for (i=0;i<dimension*2;i+=scalar)
	{
		if (i%2 == 0)printf("Real: ");
		else printf("Imag: ");
		for (j=0;j<dimension*2;j+=2) printf("%5.2f|",(hamiltonian[j*dimension+i]+0.0));
		printf("\n");
	}
	printf("\n");
}
