#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex>

#define NEIGHBORS 4
#define constant 2
#define CANT acos(1.0/3.0)/2.0
#define EIGENVALUES 40
#define nx 3
#define ny 3
#define nu 2/(nx*ny)
#define V 1.1
#define T 4.1
#define K 3
#define J constant*2
#define B 1
#define dT 0.5
#define Tao 10
#define N Tao/dT
/*todo:
  evolution function
  get matrix in float form
  evolution function
    -diagonalize
    -functionalize
  cost function
    -complex conjugate, should I have complex numbers?
    -hYdotY?
  documentation
*/

using namespace std;

typedef complex <double> dcmplx;
typedef struct
{
    int coordX, coordY;
} Coordinates;
typedef struct
{
   int i, j;
   dcmplx v;
} Telt;
const dcmplx I(0,1);

float BList[][nx*ny]=
{{0.0,1.1,2.2,3.3,4.4,5.5,6.6,7.7,8.8},{1.1,2.2,3.3,4.4,5.5,6.6,7.7,8.8,9.9},{2.2,3.3,4.4,5.5,6.6,7.7,8.8,9.9,11.0},{3.3,4.4,5.5,6.6,7.7,8.8,9.9,11.0,12.1},
{4.4,5.5,6.6,7.7,8.8,9.9,11.0,12.1,13.2},{5.5,6.6,7.7,8.8,9.9,11.0,12.1,13.2,14.3},{6.6,7.7,8.8,9.9,11.0,12.1,13.2,14.3,15.4},{7.7,8.8,9.9,11.0,12.1,13.2,14.3,15.4,16.5},
{8.8,9.9,11.0,12.1,13.2,14.3,15.4,16.5,17.6},{9.9,11.0,12.1,13.2,14.3,15.4,16.5,17.6,18.7},{11.0,12.1,13.2,14.3,15.4,16.5,17.6,18.7,19.8},{12.1,13.2,14.3,15.4,16.5,17.6,18.7,19.8,20.9},
{13.2,14.3,15.4,16.5,17.6,18.7,19.8,20.9,22.0},{14.3,15.4,16.5,17.6,18.7,19.8,20.9,22.0,23.1},{15.4,16.5,17.6,18.7,19.8,20.9,22.0,23.1,24.2},{16.5,17.6,18.7,19.8,20.9,22.0,23.1,24.2,25.3},
{17.6,18.7,19.8,20.9,22.0,23.1,24.2,25.3,26.4},{18.7,19.8,20.9,22.0,23.1,24.2,25.3,26.4,27.5},{19.8,20.9,22.0,23.1,24.2,25.3,26.4,27.5,28.6},{20.9,22.0,23.1,24.2,25.3,26.4,27.5,28.6,29.7}};

float KList[][nx*ny*2]=
{{0.0,1.1,2.2,3.3,4.4,5.5,6.6,7.7,8.8,9.9,11.0,12.1,13.2,14.3,15.4,16.5,17.6,18.7},{1.1,2.2,3.3,4.4,5.5,6.6,7.7,8.8,9.9,11.0,12.1,13.2,14.3,15.4,16.5,17.6,18.7,19.8},{2.2,3.3,4.4,5.5,6.6,7.7,8.8,9.9,11.0,12.1,13.2,14.3,15.4,16.5,17.6,18.7,19.8,20.9},{3.3,4.4,5.5,6.6,7.7,8.8,9.9,11.0,12.1,13.2,14.3,15.4,16.5,17.6,18.7,19.8,20.9,22.0},
{4.4,5.5,6.6,7.7,8.8,9.9,11.0,12.1,13.2,14.3,15.4,16.5,17.6,18.7,19.8,20.9,22.0,23.1},{5.5,6.6,7.7,8.8,9.9,11.0,12.1,13.2,14.3,15.4,16.5,17.6,18.7,19.8,20.9,22.0,23.1,24.2},{6.6,7.7,8.8,9.9,11.0,12.1,13.2,14.3,15.4,16.5,17.6,18.7,19.8,20.9,22.0,23.1,24.2,25.3},{7.7,8.8,9.9,11.0,12.1,13.2,14.3,15.4,16.5,17.6,18.7,19.8,20.9,22.0,23.1,24.2,25.3,26.4},
{8.8,9.9,11.0,12.1,13.2,14.3,15.4,16.5,17.6,18.7,19.8,20.9,22.0,23.1,24.2,25.3,26.4,27.5},{9.9,11.0,12.1,13.2,14.3,15.4,16.5,17.6,18.7,19.8,20.9,22.0,23.1,24.2,25.3,26.4,27.5,28.6},{11.0,12.1,13.2,14.3,15.4,16.5,17.6,18.7,19.8,20.9,22.0,23.1,24.2,25.3,26.4,27.5,28.6,29.7},{12.1,13.2,14.3,15.4,16.5,17.6,18.7,19.8,20.9,22.0,23.1,24.2,25.3,26.4,27.5,28.6,29.7,30.8},
{13.2,14.3,15.4,16.5,17.6,18.7,19.8,20.9,22.0,23.1,24.2,25.3,26.4,27.5,28.6,29.7,30.8,31.9},{14.3,15.4,16.5,17.6,18.7,19.8,20.9,22.0,23.1,24.2,25.3,26.4,27.5,28.6,29.7,30.8,31.9,33.0},{15.4,16.5,17.6,18.7,19.8,20.9,22.0,23.1,24.2,25.3,26.4,27.5,28.6,29.7,30.8,31.9,33.0,34.1},{16.5,17.6,18.7,19.8,20.9,22.0,23.1,24.2,25.3,26.4,27.5,28.6,29.7,30.8,31.9,33.0,34.1,35.2},
{17.6,18.7,19.8,20.9,22.0,23.1,24.2,25.3,26.4,27.5,28.6,29.7,30.8,31.9,33.0,34.1,35.2,36.3},{18.7,19.8,20.9,22.0,23.1,24.2,25.3,26.4,27.5,28.6,29.7,30.8,31.9,33.0,34.1,35.2,36.3,37.4},{19.8,20.9,22.0,23.1,24.2,25.3,26.4,27.5,28.6,29.7,30.8,31.9,33.0,34.1,35.2,36.3,37.4,38.5},{20.9,22.0,23.1,24.2,25.3,26.4,27.5,28.6,29.7,30.8,31.9,33.0,34.1,35.2,36.3,37.4,38.5,39.6}};


float JList[][nx*ny*2]=
{{7.0,8.1,9.2,10.3,11.4,12.5,13.6,14.7,15.8,16.9,18.0,19.1,20.2,21.3,22.4,23.5,24.6,25.7},{8.1,9.2,10.3,11.4,12.5,13.6,14.7,15.8,16.9,18.0,19.1,20.2,21.3,22.4,23.5,24.6,25.7,26.8},{9.2,10.3,11.4,12.5,13.6,14.7,15.8,16.9,18.0,19.1,20.2,21.3,22.4,23.5,24.6,25.7,26.8,27.9},{10.3,11.4,12.5,13.6,14.7,15.8,16.9,18.0,19.1,20.2,21.3,22.4,23.5,24.6,25.7,26.8,27.9,29.0},
{11.4,12.5,13.6,14.7,15.8,16.9,18.0,19.1,20.2,21.3,22.4,23.5,24.6,25.7,26.8,27.9,29.0,30.1},{12.5,13.6,14.7,15.8,16.9,18.0,19.1,20.2,21.3,22.4,23.5,24.6,25.7,26.8,27.9,29.0,30.1,31.2},{13.6,14.7,15.8,16.9,18.0,19.1,20.2,21.3,22.4,23.5,24.6,25.7,26.8,27.9,29.0,30.1,31.2,32.3},{14.7,15.8,16.9,18.0,19.1,20.2,21.3,22.4,23.5,24.6,25.7,26.8,27.9,29.0,30.1,31.2,32.3,33.4},
{15.8,16.9,18.0,19.1,20.2,21.3,22.4,23.5,24.6,25.7,26.8,27.9,29.0,30.1,31.2,32.3,33.4,34.5},{16.9,18.0,19.1,20.2,21.3,22.4,23.5,24.6,25.7,26.8,27.9,29.0,30.1,31.2,32.3,33.4,34.5,35.6},{18.0,19.1,20.2,21.3,22.4,23.5,24.6,25.7,26.8,27.9,29.0,30.1,31.2,32.3,33.4,34.5,35.6,36.7},{19.1,20.2,21.3,22.4,23.5,24.6,25.7,26.8,27.9,29.0,30.1,31.2,32.3,33.4,34.5,35.6,36.7,37.8},
{20.2,21.3,22.4,23.5,24.6,25.7,26.8,27.9,29.0,30.1,31.2,32.3,33.4,34.5,35.6,36.7,37.8,38.9},{21.3,22.4,23.5,24.6,25.7,26.8,27.9,29.0,30.1,31.2,32.3,33.4,34.5,35.6,36.7,37.8,38.9,40.0},{22.4,23.5,24.6,25.7,26.8,27.9,29.0,30.1,31.2,32.3,33.4,34.5,35.6,36.7,37.8,38.9,40.0,41.1},{23.5,24.6,25.7,26.8,27.9,29.0,30.1,31.2,32.3,33.4,34.5,35.6,36.7,37.8,38.9,40.0,41.1,42.2},
{24.6,25.7,26.8,27.9,29.0,30.1,31.2,32.3,33.4,34.5,35.6,36.7,37.8,38.9,40.0,41.1,42.2,43.3},{25.7,26.8,27.9,29.0,30.1,31.2,32.3,33.4,34.5,35.6,36.7,37.8,38.9,40.0,41.1,42.2,43.3,44.4},{26.8,27.9,29.0,30.1,31.2,32.3,33.4,34.5,35.6,36.7,37.8,38.9,40.0,41.1,42.2,43.3,44.4,45.5},{27.9,29.0,30.1,31.2,32.3,33.4,34.5,35.6,36.7,37.8,38.9,40.0,41.1,42.2,43.3,44.4,45.5,46.6}};


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



void construct_device_hamiltonian(int *table, unsigned long long int *b,int electrons,int dimension, int *ham_dev, int sites, int lattice[][nx], int index);
void construct_model_hamiltonian(int *table, unsigned long long int *b,int electrons,int dimension, int *ham_mod, int sites,  int lattice[][nx]);
int find(int dimension,int k,unsigned long long int *v, unsigned long long int *b);
unsigned long long int choose(int n, int k);
int combinations ( int n, int k, unsigned long long int *b,int *tab);
void construct_lattice(int lattice[][nx]);
Coordinates find_coordinates(int site_num, int lattice[][nx]);
int hop(int k, unsigned long long int *v1, unsigned long long int *v2,int n, int j);
void calc_sigma_zK(unsigned long long int b, int *ham_mod, int lattice[][nx], int dimension, int i, int index);
void calc_sigma_zV(unsigned long long int state, int *hamiltonian, int lattice[][nx], int dimension, int i);
void calc_B(unsigned long long int state, int *ham_dev, int dimension, int i);
void print_hamiltonian(int* hamiltonian, int dimension);
float evolve(int *table, unsigned long long int *b,int electrons,int dimension, int sites,  int lattice[][nx], int *ham_dev);
float functionalize(int state, int *ham);
float cost(int state, int *ham_mod);
void printArray(int scalar, int offset);



int main (int argc, char *argv[])
{
   int *table,*v, *v2,i,j,k,nev,sites, electrons, dimension,z,x,*ham_mod,*ham_dev;
   unsigned long long int *b;
   float state;
   int lattice[nx][nx];
   int expectation = 0;

   sites=nx*nx;
   electrons=sites*nu; // 1
   dimension=choose(sites,electrons);
   //printArray(1,0);
   //printArray(2,0);
   //printArray(2,7);

   b=new unsigned long long int[dimension];
   table=(int*) malloc(electrons*dimension*sizeof(int));
   ham_mod = (int *) malloc (dimension*dimension*sizeof (int));
   ham_dev = (int *) malloc (dimension*dimension*sizeof (int));
   ham_mod[0]=0;
   ham_mod[1]=0;
   ham_dev[0]=0;
   ham_dev[1]=0;

   construct_lattice(lattice);
   combinations (sites,electrons,b,table);
   construct_model_hamiltonian(table, b, electrons, dimension, ham_mod, sites, lattice);
   state = evolve(table,b,electrons,dimension,sites,lattice, ham_dev);
   expectation = cost(state, ham_mod);
   //printf("%i\n",expectation);
   //printf("The Model Hamiltonian:\n");
   //print_hamiltonian(ham_mod, dimension);
   //printf("The Device Hamiltonian:\n");
   print_hamiltonian(ham_dev, dimension);
   exit (0);
}


float evolve(int *table, unsigned long long int *b,int electrons,int dimension, int sites,  int lattice[][nx], int *ham_dev)
{
  int i;
  float state;
  state = 0; //starting state

  for (i=0; i<N;i++)
  {
    memset(ham_dev, 0, sizeof(ham_dev)); //This should get replaced when converting to flaot
    construct_device_hamiltonian(table, b, electrons, dimension, ham_dev, sites,lattice, i);
    //printf("state before mult :%i\n", state); //checking the state
    state = functionalize(state, ham_dev);
    //printf("state AFTER mult :%i\n\n", state); //checking the state
  }
  return state;
}


float functionalize(int state, int *ham)
{
    //e^(i*dt*ham)*state
    //return state;
    return 0;
}


float cost(int state, int *ham_mod)
{
  //(ham_mod*state)dot(state)
  return 0;
}


void construct_device_hamiltonian(int *table, unsigned long long int *b,int electrons,int dimension, int *ham_dev, int sites, int lattice[][nx], int index)
/*Constructing the hamiltonian matrix for the device hamiltonian

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
sites : integer
      total number of sites in the lattice (nx^2)

Updates
-------
Tdev : array-pointer
    Stores the values of the model-hamiltonian matrix*/
{
  int i,state,j, k,l,n,x,y,z,p ,site,neighbor,sign,bond;
  unsigned long long int *v1,*v2;
  Coordinates coord;
  v1=(unsigned long long int*) malloc(sizeof(unsigned long long int));
  v2=(unsigned long long int*) malloc(sizeof(unsigned long long int));

  for (i=0;i<dimension;i++)
  {
    for (j=0;j<electrons;j++)
    {
      site=table[i*electrons+j];
      coord = find_coordinates(site, lattice);
      x = coord.coordX;
      y = coord.coordY;
      int neighbors[NEIGHBORS] = {lattice[(x+1)%nx][y], lattice[(x+(nx-1))%nx][y], lattice[x][(y+1)%nx], lattice[x][(y+(nx-1))%nx]};

      for (z=0; z<NEIGHBORS; z++)
      {

        if (((1ULL<<(neighbors[z]-1))&b[i])==0)//making sure neighbor is not occupied, otherwise nothing happens
        {

          memcpy(&v1,&b[i], sizeof(unsigned long long int));
          hop(electrons, v1, v2,site, neighbors[z]);
          state=find(dimension,electrons,v2, b);
          bond = find_bond[site-1][neighbors[z]-1];
          if (bond<=0 || bond>18) printf("something bad is happening, bond = %i\n", bond);
          //printf("pre = %.1f\n", JList[index][bond]);
          ham_dev[dimension*i+state-1] += JList[index][bond];
          //printf("post = %.1f\n\n", ham_dev[dimension*i+state-1]);

        }
      }
    }
    calc_sigma_zK(b[i], ham_dev,lattice,dimension, i, index);//The K term calculation
    for (p=0; p<nx*ny;p++)//The B term calculation
    {
      if((1ULL<<p)&b[i]) ham_dev[(dimension*i)+i]+=BList[index][p];
      else ham_dev[(dimension*i)+i]-=BList[index][p];
    }
  }
}




void construct_model_hamiltonian(int *table, unsigned long long int *b,int electrons,int dimension, int *ham_mod, int sites,  int lattice[][nx])
/*Constructing the hamiltonian matrix for the model hamiltonian

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
sites : integer
      total number of sites in the lattice (nx^2)

Updates
-------
Tmod : array-pointer
    Stores the values of the model-hamiltonian matrix*/
{
  int i,state,j, k,l,n,x,y,z ,site,neighbor,sign;
  unsigned long long int *v1,*v2;
  Coordinates coord;
  v1=(unsigned long long int*) malloc(sizeof(unsigned long long int));
  v2=(unsigned long long int*) malloc(sizeof(unsigned long long int));

  for (i=0;i<dimension;i++)
  {

    for (j=0;j<electrons;j++)
    {
      site=table[i*electrons+j];
      coord = find_coordinates(site, lattice);
      x = coord.coordX;
      y = coord.coordY;
      int neighbors[NEIGHBORS] = {lattice[(x+1)%nx][y], lattice[(x+(nx-1))%nx][y], lattice[x][(y+1)%nx], lattice[x][(y+(nx-1))%nx]};

      for (z=0; z<NEIGHBORS; z++)
      {

        if (((1ULL<<(neighbors[z]-1))&b[i])==0)//checking if the neighbor is occupied
        {

          memcpy(&v1,&b[i], sizeof(unsigned long long int));
          sign = hop(electrons, v1, v2,site, neighbors[z]);
          if (sign==0) sign=1;
          else sign=-1;
          state=find(dimension,electrons,v2, b);
          ham_mod[dimension*i+state-1] -= (T*sign);
        }
      }
    }
    calc_sigma_zV(b[i], ham_mod,lattice,dimension, i);
  }
}


void calc_sigma_zK(unsigned long long int state, int *ham_dev, int lattice[][nx], int dimension, int i, int index)
{
  int p,n,sign,x,y,site,bond;
  unsigned long long int comparison;
  Coordinates coord;
  for (n=1;n<(nx*ny);n++)//Only need nx*ny-1 comparisons, the last site has already had all neighbors compared
  {
    site=n;
    coord = find_coordinates(site, lattice);
    x = coord.coordX;
    y = coord.coordY;
    int neighbors[NEIGHBORS] = {lattice[(x+1)%nx][y], lattice[(x+(nx-1))%nx][y], lattice[x][(y+1)%nx], lattice[x][(y+(nx-1))%nx]};

    for (p=0; p<4;p++)
    {
      if (neighbors[p] > site)
      {
        sign = -1;
        comparison = (1ULL<<(neighbors[p]-1))+(1ULL<<(site-1));
        if((comparison&state)==comparison || (comparison&state)==0) sign = 1;
        bond = find_bond[site][neighbors[p]];
        ham_dev[(dimension*i)+i] += sign*(KList[index][bond]);
      }
    }
  }
}



void calc_sigma_zV(unsigned long long int state, int *ham_mod, int lattice[][nx], int dimension, int i)
{
  int p,n,sign,x,y,site,term;
  unsigned long long int comparison;
  Coordinates coord;

  for (n=1;n<(nx*ny);n++)//Only need nx*ny-1 comparisons, the last site has already had all neighbors compared
  {
    site=n;
    coord = find_coordinates(site, lattice);
    x = coord.coordX;
    y = coord.coordY;
    int neighbors[NEIGHBORS] = {lattice[(x+1)%nx][y], lattice[(x+(nx-1))%nx][y], lattice[x][(y+1)%nx], lattice[x][(y+(nx-1))%nx]};

    for (p=0; p<4;p++)
    {
      if (neighbors[p] > site)
      {
        sign = -1;
        comparison = (1ULL<<(neighbors[p]-1))+(1ULL<<(site-1));
        if((comparison&state)==comparison || (comparison&state)==0) sign = 1;
        ham_mod[(dimension*i)+i] += sign*V;
      }
    }
  }
}


int hop(int k, unsigned long long int *v1, unsigned long long int *v2,int n, int j)
/*given a state v1 (a combination corresponding to k occupied sites) generates the state
v2 obtained by hopping from site v1[n] to site j (which should not be in v1), and outputs the fermionic sign

Parameters
----------
k : integer
    The total number of electrons per state
v1 : array-pointer (vector)
    A state, a combination corresponding to k occupied sites
v2 : array-pointer (vector)
    The new resulting state after the hop
n : integer
    the hop from site
j : integer
    the hop-to site

Updates
-------
v2 : array-pointer (vector)
    The new resulting state after the hop

Returns
-------
i%2 : integer
    The fermionic sign from the resulting hop*/
{
    unsigned long long int p,i,a,vtemp;
    int z_count = 0;
    vtemp = (uintptr_t) v1;
    p = (1ULL << (n-1)) + (1ULL << (j-1));
    for (i=n;i<j-1;i++)  if((1ULL<<i) & (vtemp)) z_count++;
    if((1ULL<<4) & (1ULL<<4));
    a = (p ^ vtemp);
    memcpy(v2, &a, sizeof(unsigned long long int));
    return z_count%2;
}

int find(int dimension,int k, unsigned long long int *v,unsigned long long int *b)
/*find the position of a given combination v (vector of length k) in the table of all combinations tab (array of length C(n,k)*k)

Parameters
----------
dimension : integer
    total number of states, n-choose-k where n is total # of electrons, and k is total number of sites.
k : integer
    number of electrons
v : array-pointer (vector)
    the state (combination of occupied sites) that we're looking for in the table
tab : array-pointer
    An array containing the location of the occupied sites. Ex: if the first two states = 110000000 and 000000011, the first 4
    entries of the table read: 1, 2, 8, 9.k : integer
b : array-pointer
    An array containing binary numbers, representing occupied (1) and unoccupied (0) sites. Ex: if 2 of 9 sites are occupied,
    sites 2 and 4, this would be repesented by 010100000. b contains a total of n-choose-k elements.

Returns
-------
mid : integer
    the location of the state in the table*/
{
   int i, first, last, mid;
   unsigned long long int vtemp;
   vtemp = *v;
   first=0;
   last=dimension-1;

   while (first <= last)
   {
      mid = (int) ((first + last) / 2.0);
      if (vtemp > b[mid]) first = mid + 1;
      else if (vtemp < b[mid]) last = mid - 1;
      else return mid+1;
   }
}



unsigned long long int choose(int n, int k)
/*Returning the number of combinations of k out of n. (n-choose-k), finding number of unordered combinations with n electrons in k sites.

Parameters
----------
n : integer
    The total number of electrons, the number being chosen.
k : integer
    The total number of sites, the number of possible choices

Returns
-------
c : integer
    The number of ways n electron can be arranged in k sites*/
{
   int i;
   unsigned long long int c;
   c=1ULL;
   for (i=0;i<k;i++) c=c*(n-i);
   for (i=0;i<k;i++) c=c/(i+1);
   return c;
}



int combinations ( int n, int k, unsigned long long int *b,int *tab)
/*Returning the number of combinations of k out of n. (n-choose-k), finding number of unordered combinations with n electrons in k sites.

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
    entries of the table read: 1, 2, 8, 9.*/
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

void construct_lattice(int lattice[][nx])
/*Constructing the lattice of dimension nx*nx

Parameters
----------
lattice[][nx] : array
    an empty array containing each site, which gets filled with a site number*/
{
  int x,i;
  for (x=0;x<nx;x++)
  {
    if (x%2 ==1)
    {
      for (i=0;i<nx;i++) lattice[x][i] = (nx*x)+i+1;
    }
    else
    {
      for (i=0;i<nx;i++) lattice[x][nx-1-i] = nx*x+i+1;
    }
  }
}

Coordinates find_coordinates(int site_num, int lattice[][nx])
/*Finding the coordinates of a given site

Parameters
----------
site_num : integer
    the site number of the site we're looking for
lattice[][nx] : array
    an array containing each site

Returns
-------
coords : Coordinates
    the x and y coordinates of the site on the lattice*/
{
  int x,i;
  for (i=0;i<nx;i++)
  {
    for (x=0;x<nx;x++)
    {

      if(site_num == lattice[i][x])
      {
        Coordinates coords = { i,x };
        return coords;
      }
    }
  }
}


void print_hamiltonian(int* hamiltonian, int dimension)
{
  int i,j;

  for (i=0;i<dimension;i++){
    for (j=0;j<dimension;j++) printf("%i ",(hamiltonian[j*dimension+i]));
    printf("\n");

  }
  printf("\n");

}

void printArray(int scalar, int offset)
{
  int z,x;
  printf("{");
  for (z=0;z<N;z++)
  {
    printf("{%.1f",z*1.1+offset);
    for (x=1;x<nx*ny*scalar;x++)
    {
      printf(",%.1f",((z+x)*1.1+offset));

    }
    if (z<N-1)
    {
      if ((z%4)==3)
      {
        printf("},\n");
      }
      else printf("},");
    }
    else printf("}};\n\n");
  }

}
