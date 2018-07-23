#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex>

#define COORD 14
#define CANT acos(1.0/3.0)/2.0
#define INTERACTION 0.35
#define nu 1/9
#define DISCRET 0
#define EIGENVALUES 40
#define T1 -0.13
#define T2 -0.02
#define T3 0.05

/*todo:
  H-device function
  H-model function
  Add hop function
      What is j and n?
*/


using namespace std;

typedef complex <double> dcmplx;
typedef struct
{
     int neighbors[COORD];
     dcmplx t[COORD];
     double V[COORD];
} Site;
typedef struct
{
   int i, j;
   dcmplx v;
} Telt;
const dcmplx I(0,1);




int construct_device_hamiltonian(int *table, unsigned long long int *b,int electrons,int dimension, Telt *Tdev);
int construct_model_hamiltonian(int *table, unsigned long long int *b,int electrons,int dimension, Telt *Tmod);
int hop(int k, int *v1, int *v2,int n, int j);
int find(int dimension,int k, int *v, int *T, unsigned long long int *b);
unsigned long long int choose(int n, int k);
int combinations ( int n, int k, unsigned long long int *b,int *tab);






int main (int argc, char *argv[])
{
   int *table,*v, *v2,i,j,k,sign, nev,tlen,Nx,Ny,sites, electrons, dimension;
   unsigned long long int *b;
   FILE *fp;
   Telt *Tdev, *Tmod;
   dcmplx **Evecs, *Evals;
   Nx = atoi(argv[1]);
   Ny = atoi(argv[2]);
   sites=Nx*Ny;
   //sites = 9;
   electrons=sites*nu; // 2
   dimension=choose(sites,electrons);

   b=new unsigned long long int[dimension];
   Tdev = (Telt *) malloc (dimension*(electrons*COORD+1)*sizeof (Telt));
   Tmod = (Telt *) malloc (dimension*(electrons*COORD+1)*sizeof (Telt));
   table=(int*) malloc(electrons*dimension*sizeof(int));
   v=(int*) malloc(electrons*sizeof(int));
   v2=(int*) malloc(electrons*sizeof(int));

/*
   nev=EIGENVALUES;
   Evals = new dcmplx[nev];
   Evecs = new dcmplx*[dimension];
   for (i=0; i<dimension; i++) Evecs[i] = new dcmplx[nev];
*/
   combinations (sites,electrons,b,table);
   construct_device_hamiltonian(table, b, electrons, dimension, Tdev);
   construct_model_hamiltonian(table, b, electrons, dimension, Tmod);
   ///tlen=construct_hamiltonian(table,b, electrons, dimension, T);
   printf("Missions success: dimnesion = %i", dimension);
   exit (0);
}




int construct_device_hamiltonian(int *table, unsigned long long int *b,int electrons,int dimension, Telt *Tdev)
/*
Constructing the hamiltonian matrix for the device hamiltonian

Parameters
----------
table : array-pointer
    An array containing the location of the occupied sites. Ex: if the first two states = 110000000 and 000000011, the first 4
    entries of the table read: 1, 2, 8, 9.k : integer
b : array-pointer
    An array containing binary numbers, representing occupied (1) and unoccupied (0) sites. Ex: if 2 of 9 sites are occupied,
    sites 2 and 4, this would be repesented by 010100000. b contains a total of n-choose-k elements.
electron : integer
      number of electrons
dimension : integer
      total number of states, n-choose-k where n is total # of electrons, and k is total number of sites.
Tdev : array-pointer (empty)

Updates
-------
Tdev : array-pointer
    Stores the values of the model-hamiltonian matrix
*/
{
  //sign=hop(electrons, v1, v2,j, neighbor+1);
  //state=find(dimension,electrons,v2, table, b);


}




int construct_model_hamiltonian(int *table, unsigned long long int *b,int electrons,int dimension, Telt *Tmod)
/*
Constructing the hamiltonian matrix for the model hamiltonian

Parameters
----------
table : array-pointer
    An array containing the location of the occupied sites. Ex: if the first two states = 110000000 and 000000011, the first 4
    entries of the table read: 1, 2, 8, 9.k : integer
b : array-pointer
    An array containing binary numbers, representing occupied (1) and unoccupied (0) sites. Ex: if 2 of 9 sites are occupied,
    sites 2 and 4, this would be repesented by 010100000. b contains a total of n-choose-k elements.
electron : integer
      number of electrons
dimension : integer
      total number of states, n-choose-k where n is total # of electrons, and k is total number of sites.
Tmod : array-pointer (empty)

Updates
-------
Tmod : array-pointer
    Stores the values of the model-hamiltonian matrix
*/
{

}



int hop(int k, int *v1, int *v2,int n, int j)
/*
given a state v1 (a combination corresponding to k occupied sites) generates the state v2 obtained by hopping from site v1[n] to site j (which should not be in v1), and outputs the fermionic sign

Parameters
----------
k : integer
    The total number of electrons per state
v1 : array-pointer (vector)
    A state, a combination corresponding to k occupied sites
v2 : array-pointer (vector)
    The new resulting state after the hop
n : integer
    TO BE DETERMINED
j : integer
    TO BE DETERMINED

Updates
-------
v2 : array-pointer (vector)
    The new resulting state after the hop

Returns
-------
i%2 : integer
    The fermionic sign from the resulting hop
*/
{
   int i,d;

   i=0;
   if (j>v1[n]) {while((n+i+1<k)&&(j>v1[n+i+1]))i++;}

   if (j<v1[n]) {while((n+i-1>=0)&&(j<v1[n+i-1]))i--;}
   memcpy(&v1[n],&v1[n+1], (k-n-1)*sizeof(int));
   memcpy(v2,v1, (n+i)*sizeof(int));
   v2[n+i]=j;
   memcpy(&v2[n+i+1],&v1[n+i], (k-n-i-1)*sizeof(int));
   return i%2;
}




int find(int dimension,int k, int *v, int *T, unsigned long long int *b)
/*
find the position of a given combination v (vector of length k) in the table of all combinations T (array of length C(n,k)*k)

Parameters
----------
dimension : integer
    total number of states, n-choose-k where n is total # of electrons, and k is total number of sites.
k : integer
    number of electrons
v : array-pointer (vector)
    the state (combination of occupied sites) that we're looking for in the table
T : array-pointer
    An array containing the location of the occupied sites. Ex: if the first two states = 110000000 and 000000011, the first 4
    entries of the table read: 1, 2, 8, 9.k : integer
b : array-pointer
    An array containing binary numbers, representing occupied (1) and unoccupied (0) sites. Ex: if 2 of 9 sites are occupied,
    sites 2 and 4, this would be repesented by 010100000. b contains a total of n-choose-k elements.

Returns
-------
mid : integer
    the location of the state in the table
*/
{
   int i, first, last, mid;
   unsigned long long int element;
   element=0ULL;

   for (i=0;i<k;i++) element=element+(1ULL << (v[i]-1));
   //for (i=0;i<k;i++) cout<<v[i]<<" ";
   //cout<<element<<"\n";
   first=0;
   last=dimension-1;


   while (first <= last)
   {
           mid = (int) ((first + last) / 2.0);
           if (element > b[mid])
               first = mid + 1;
           else if (element < b[mid])
               last = mid - 1;
              else
              {
        return mid;}

   }
}



unsigned long long int choose(int n, int k)
/*
Returning the number of combinations of k out of n. (n-choose-k), finding number of unordered combinations with n electrons in k sites.

Parameters
----------
n : integer
    The total number of electrons, the number being chosen.
k : integer
    The total number of sites, the number of possible choices

Returns
-------
c : integer
    The number of ways n electron can be arranged in k sites
*/
{
   int i;
   unsigned long long int c;
   c=1ULL;
   for (i=0;i<k;i++) c=c*(n-i);
   for (i=0;i<k;i++) c=c/(i+1);
   return c;
}



int combinations ( int n, int k, unsigned long long int *b,int *tab)
/*
Returning the number of combinations of k out of n. (n-choose-k), finding number of unordered combinations with n electrons in k sites.

Parameters
----------
n : integer
    The total number of electrons.
k : integer
    The total number of sites.
b : array-pointer (empty)
tab : array-pointer (empty)

Updates
-------
b : array
    An array containing binary numbers, representing occupied (1) and unoccupied (0) sites. Ex: if 2 of 9 sites are occupied,
    sites 2 and 4, this would be repesented by 010100000. b contains a total of n-choose-k elements.
tab : array
    An array containing the location of the occupied sites. Ex: if the first two states = 110000000 and 000000011, the first 4
    entries of the table read: 1, 2, 8, 9.
*/
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
