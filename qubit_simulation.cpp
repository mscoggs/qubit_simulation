#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <complex>
#include <math.h>

#define NEIGHBORS 4
#define nx 3
#define ny 3
#define nu 2/(nx*ny)
#define V 1.1
#define T 4.1
#define JJ 2.6
#define KK 2.7
#define BB 2.8
#define dT 0.5
#define Tao 10
#define PsiStartState 2
#define Ncount Tao/dT

/*todo:
  evolution function
    -exp needs to work w/ imaginary
  cost function

  documentation
  algorithmic way to find bonds
*/


using namespace std;
typedef struct
{
    int coordX, coordY;
} Coordinates;


int find_bond[][nx*ny]=
{{0,1,3,0,0,10,16,0,0},
{1,0,2,0,11,0,0,17,0},
{3,2,0,12,0,0,0,0,18},
{0,0,12,0,5,6,0,0,15},
{0,11,0,5,0,4,0,14,0},
{10,0,0,6,4,0,13,0,0},
{16,0,0,0,0,13,0,7,9},
{0,17,0,0,14,0,7,0,8},
{0,0,18,15,0,0,9,8,0}};



double evolve(int *table, unsigned long long int *b,int electrons,int dimension, int sites,  int lattice[][nx], double *ham_dev);
double cost(int state, double *ham_mod);
void construct_device_hamiltonian(int *table, unsigned long long int *b,int electrons,int dimension, double *ham_dev, int sites, int lattice[][nx],double *k_list, double *j_list, double *b_list);
void construct_model_hamiltonian(int *table, unsigned long long int *b,int electrons,int dimension, double *ham_mod, int sites,  int lattice[][nx]);
void diag_hermitian_real_double(int N,  double *A, double *Vdag,double *D);
void exp_general_complex_double(int N, double *A, double *B);
int hop(unsigned long long int *v1, unsigned long long int *v2,int n, int j);
int find(int dimension,unsigned long long int *v, unsigned long long int *b);
unsigned long long int choose(int n, int k);
int combinations ( int n, int k, unsigned long long int *b,int *tab);
void init_lists(double *k_list, double *j_list,double *b_list);
void change_lists(double *k_list, double *j_list,double *b_list, int index);
void construct_lattice(int lattice[][nx]);
Coordinates find_coordinates(int site_num, int lattice[][nx]);
void print_hamiltonian(double* hamiltonian, int dimension, char print_imaginary);
extern "C" int dsyev_(char *JOBZ, char *UPLO, int *N, double *Vdag, int *LDA, double *D, double *WORK, int *LWORK, int *INFO);
void print_hamiltonianR(double *hamiltonian, int dimension);
void exp_diaganolized_mat(double *ham_real, double *Vdag, double* D, int dimension);





int main (int argc, char *argv[])
{
   int *table,*v,*v2,expectation,sites,electrons,dimension,lattice[nx][nx],*j_list,*b_list,*k_list;
   unsigned long long int *b;
   double psi,*ham_mod,*ham_dev;
   sites=nx*ny;
   electrons=sites*nu;
   dimension=choose(sites,electrons);

   b=new unsigned long long int[dimension];
   table=(int*) malloc(electrons*dimension*sizeof(int));
   ham_mod = (double *) malloc (2*dimension*dimension*sizeof (double));
   ham_dev = (double *) malloc (2*dimension*dimension*sizeof (double));

   construct_lattice(lattice);
   combinations (sites,electrons,b,table);
   construct_model_hamiltonian(table, b, electrons, dimension, ham_mod, sites, lattice);
   psi = evolve(table,b,electrons,dimension,sites,lattice, ham_dev);
   expectation = cost(psi, ham_mod);
   printf("The Device Hamiltonian:\n");
   //print_hamiltonian(ham_dev, dimension, 'n');
   //printf("The Model Hamiltonian:\n");
   //print_hamiltonian(ham_mod, dimension, 'n');
   exit (0);
}


double evolve(int *table, unsigned long long int *b,int electrons,int dimension, int sites,  int lattice[][nx], double *ham_dev)

{
  int i,j;
  double *psi1,*psi2,*ham_t_i, *ham_diag,*exp_matrix,*D, *Vdag,*b_list,*k_list,*j_list;
  k_list = (double *) malloc(2*nx*ny*sizeof(double));
  j_list = (double *) malloc(2*nx*ny*sizeof(double));
  b_list = (double *) malloc(2*nx*ny*sizeof(double));
  ham_diag = (double *) malloc(dimension*dimension*sizeof(double));
  D = (double*) malloc(sizeof(double)*dimension);
  Vdag = (double*) malloc(sizeof(double)*dimension*dimension);
  ham_t_i = (double *) malloc (2*dimension*dimension*sizeof (double));
  exp_matrix = (double *) malloc (2*dimension*dimension*sizeof (double));
  psi1 = (double *) malloc (dimension*sizeof(double));
  psi2 = (double *) malloc (dimension*sizeof(double));

  for (i=0; i<dimension;i++) psi1[i] = 0, psi2[i] = 0;


  psi1[PsiStartState] = 1; //arbitrary starting state, convert to vector
  psi2[PsiStartState] = 1; //arbitrary starting state, convert to vector

  init_lists(k_list, j_list, b_list);

  for (i=Ncount-1; i<Ncount;i++)
  {

    change_lists(k_list,j_list,b_list,i);
    construct_device_hamiltonian(table, b, electrons, dimension, ham_dev, sites,lattice, k_list, j_list, b_list);

    for (j=0; j<dimension*dimension*2-1; j++)
      exp_matrix[j+1] = 0;
      ham_t_i[j+1] = (ham_dev[j]*-dT)+0.0; //multiplying by -i*dt
    }
    ham_t_i[0]=0;
    exp_matrix[0] = 0;

    for (j =0; j<dimension*dimension; j++)
    {
	    ham_diag[j] = ham_dev[2*j];
    }


    diag_hermitian_real_double(dimension, ham_diag,Vdag, D);
    exp_diaganolized_mat(ham_diag, Vdag, D, dimension);
    psi1 = 0; //ham_real*psi1
    //exp_general_complex_double(dimension, ham_t_i, exp_matrix) D_mat[i*dim+i=]exp(D[i]), all else 0
    psi2 = 0; //exp_matrix*psi, in vector and matrix form

  }
  return 0;
}

void exp_diaganolized_mat(double *ham_diag, double *Vdag, double* D, int dimension)
{
   int i;
   double *exp_D;
   exp_D = (double *) malloc (dimension*dimension*sizeof (double));
   for (i =0; i<dimension*dimension; i++) exp_D[i] = 0, ham_diag[i] = 0;

   for (i =0; i<dimension; i++)
   {
    exp_D[i*dimension+i] = exp(D[i]*-dt*i);/////doesnt account for EXP to the imaginry
   }
   //matrix mult hamreal= (vDAG*expD*vdag)

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


void exp_general_complex_double(int N, double *A, double *B)
{
     int M,K,ii,jj,kk,s,p,q;
     double *row_norm,*X,*Y,*Z,*E,*D;
     double norm,c;
     int * IPIV;
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
     int LDA=N;
     int LDB=N;
     int LDC=N;
     int NRHS=N;
     int INFO;

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
         //zgemm_ (&TRANSA, &TRANSB, &M, &N, &K, ALPHA, Z, &LDA, X, &LDB, BETA, Y, &LDC); //matrix mult
         memcpy(X,Y,2*N*N*sizeof(double));

         for(ii=0;ii<2*N*N;ii++) Y[ii]=c*X[ii];
         for(ii=0;ii<2*N*N;ii++) E[ii]=E[ii]+Y[ii];
         if (p==1)for(ii=0;ii<2*N*N;ii++) D[ii]=D[ii]+Y[ii];
         else for(ii=0;ii<2*N*N;ii++) D[ii]=D[ii]-Y[ii];
         p=abs(1-p);
     }

     //zgetrf_ (&M, &N, D, &LDA, IPIV, &INFO);http://www.netlib.org/lapack/explore-3.1.1-html/zgetrf.f.html
     if (INFO !=0) printf("DIAGONALIZATION ERROR, INFO = %i\n", INFO);
     //zgetrs_( &TRANS, &N, &NRHS,D,    &LDA, IPIV,E, &LDB, &INFO);//solves system of equations http://www.netlib.org/lapack/explore-3.1.1-html/zgetrs.f.html
     if (INFO !=0) printf("DIAGONALIZATION ERROR, INFO = %i\n", INFO);

     for (kk=1;kk<=s;kk++)
     {
         memcpy(X,E,2*N*N*sizeof(double));

         //zgemm_ (&TRANSA, &TRANSB, &M, &N, &K, ALPHA, E, &LDA, X, &LDB, BETA, Y, &LDC);//matrixmultiplication
         if (INFO !=0) printf("DIAGONALIZATION ERROR, INFO = %i\n", INFO);
         memcpy(E,Y,2*N*N*sizeof(double));

     }
     memcpy(B,E,2*N*N*sizeof(double));
     free(row_norm);
     free(X);
     free(Y);
     free(E);
     free(D);;
     free(IPIV);
}






double cost(int state, double *ham_mod)
{
  //(ham_mod*state)dot(state)
  return 0;
}


void construct_device_hamiltonian(int *table, unsigned long long int *b,int electrons,int dimension, double *ham_dev, int sites, int lattice[][nx],double *k_list, double *j_list, double *b_list)
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
          ham_dev[(dimension*i+state-1)*2] += j_list[bond];
        }
      }
    }

    for (j=1;j<(nx*ny);j++)//The K term calculation
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
          ham_dev[((dimension*i)+i)*2] += k_list[bond];
        }
      }
    }

    for (j=0; j<nx*ny;j++)//The B term calculation
    {
      if((1ULL<<j)&b[i]) ham_dev[((dimension*i)+i)*2]+=b_list[j];
      else ham_dev[((dimension*i)+i)*2]-=b_list[j];
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

    for (j=1;j<(nx*ny);j++)//The V term calculation
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
   for (i=0;i<k;i++)
   {
       x=x+(1ULL<<i);
   }
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
	int i;
	for (i=0; i<nx*nx*2; i++)
	{
		j_list[i] = 0.05*i*pow(-1,i);
		k_list[i] = 0.03*i*pow(-1,i+1);
	}
	for (i=0; i<nx*nx; i++)
	{
		b_list[i] = 0.09*i*pow(-1,i);
	}
}

void change_lists(double *k_list,double *j_list,double *b_list, int index)
{
	if (index==0) return;
	int i;
	for (i=0; i<nx*nx*2; i++)
	{
		j_list[i] = j_list[i]*.9;
		k_list[i] = k_list[i]*.5;
	}
	for (i=0; i<nx*nx; i++)
	{
		b_list[i] = b_list[i]*.9*pow(-1,i);
	}



}
void construct_lattice(int lattice[][nx])
/*Constructing the lattice of dimension nx*nx*/
{
  int x,y;
  for (x=0;x<nx;x++)
  {
    if (x%2 ==1)
    {
      for (y=0;y<nx;y++) lattice[x][y] = (nx*x)+y+1;
    }
    else
    {
      for (y=0;y<nx;y++) lattice[x][nx-1-y] = nx*x+y+1;
    }
  }
}



Coordinates find_coordinates(int site_num, int lattice[][nx])
/*Finding the coordinates of a given site*/
{
  int x,y;
  for (x=0;x<nx;x++)
  {
    for (y=0;y<ny;y++)
    {

      if(site_num == lattice[x][y])
      {
        Coordinates coords = { x,y };
        return coords;
      }
    }
  }
}


void print_hamiltonianR(double* hamiltonian, int dimension)
{
  int i,j;
  for (i=0;i<dimension;i++){
    for (j=0;j<dimension;j++) printf("%.1f",(hamiltonian[j*dimension+i]));
    printf("\n");
  }
  printf("\n");
}

void print_hamiltonian(double* hamiltonian, int dimension, char print_imaginary)
{
  int i,j,scalar;
  scalar=2;
  if(print_imaginary =='y') scalar =1;
  for (i=0;i<dimension*2;i+=scalar){
    for (j=0;j<dimension*2;j+=2) printf("%.1f",(hamiltonian[j*dimension+i]));
    printf("\n");
  }
  printf("\n");
}
