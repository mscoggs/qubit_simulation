#include <iostream>
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





using namespace std;

typedef complex<double> dcmplx;
typedef struct
{
     int neighbors[COORD];
     dcmplx t[COORD];
     double V[COORD];

} Site;
typedef struct {
   int i, j;
   dcmplx v;
} Telt;
const dcmplx I(0,1);





void av(int n, dcmplx *in, dcmplx *out,int tlen, Telt *T);
#include "znaupd.h"
int construct_hamiltonian( Site* lattice, int *table, unsigned long long
int *b,int electrons,int dimension, Telt *T);
void construct_lattice( Site* lattice, int Nx, int Ny, double V,double
delta, int *ising, double phi, double tp, double tpp, double tpph);
unsigned long long int C(int n, int k);
int combinations ( int n, int k, unsigned long long int *b,int *tab) ;
int find(int dimension,int k, int *v, int *T, unsigned long long int *b);
int hop(int k, int *v1, int *v2,int n, int j);


int main (int argc, char *argv[])
{
         int     *table,*v, *v2,i,j,k,sign, nev,tlen;
     unsigned long long int *b;
     FILE *fp;
     int *ising,Nx,Ny,sites, electrons, dimension;
     Site *lattice;
     Telt *T;
     dcmplx **Evecs, *Evals;
     double phi;

         Nx = atoi (argv[1]);
     Ny = atoi (argv[2]);
     phi=atof(argv[3]);

     sites=3*Nx*Ny;
     //phi=0;
     electrons=sites*nu;
     dimension=C(sites,electrons);

     b=new unsigned long long int[dimension];
     T = (Telt *) malloc (dimension*(electrons*COORD+1)*sizeof (Telt));
     ising=(int*) malloc(sites*sizeof(int));
     table=(int*) malloc(electrons*dimension*sizeof(int));
     v=(int*) malloc(electrons*sizeof(int));
     v2=(int*) malloc(electrons*sizeof(int));
     lattice =(Site *) malloc(sites*sizeof(Site));

     nev=EIGENVALUES;
     Evals = new dcmplx[nev];
     Evecs = new dcmplx*[dimension];
     for (i=0; i<dimension; i++) Evecs[i] = new dcmplx[nev];

     for (i=0;i<sites;i++)ising[i]=1;










         combinations (sites,electrons,b,table);

     //cout <<"table done  \n";
     /*for (j=0;j<dimension;j++)
     {
         for (i=0;i<electrons;i++)cout << table[electrons*j+i]<<" ";
         cout <<" "<<b[j]<<" ";
         cout <<"\n";
     }*/

     //for (j=0; j<DISCRET; j++)
     for (j=0; j<=DISCRET; j++)
     {
     //phi=j*2.0*M_PI/DISCRET;
     phi=2.0*M_PI*phi;
     construct_lattice( lattice, Nx,Ny, INTERACTION,CANT,ising,
phi,T1,T2,T3);
     //cout <<"lattice done  \n";
     //for (j=0; j<sites; j++) {cout<<"site :"<<j<<":\n"; for (k=0;
k<COORD; k++)cout<<"neighbor["<<k<<"]="<<lattice[j].neighbors[k]<<"
t="<<lattice[j].t[k]<<"\n";}



     tlen=construct_hamiltonian( lattice,table,b, electrons, dimension, T);
     //cout <<"construction done, tlen dimension "<<tlen<<"
<<dimension<< " \n";
     //cout << "tlen dimension " << tlen<<" "<<dimension << "\n";

     /*for (j=0; j<tlen; j++)


     {    for (i=0; i<electrons; i++)
cout<<table[(T[j].i)*electrons+i]-1<<" " ;
         cout<<", ";
         for (i=0; i<electrons; i++)
cout<<table[(T[j].j)*electrons+i]-1<<" " ;
         cout<<": ";
         cout<<real(T[j].v)<<" "<<imag(T[j].v)<<"\n";
     }*/

     znaupd(dimension, nev, Evals, Evecs,tlen, T);

     cout <<phi<<" ";for (i=nev-1; i>=0; i--)  cout << real(Evals[i] )<<
" ";cout <<"\n";
     }



/*v[0]=1;
v[1]=3;
v[2]=4;
for (j=0; j<electrons; j++)std::cout << v[j]<<" ";std::cout << "\n";
sign=hop(electrons, v, v2,0, 2);
for (j=0; j<electrons; j++)std::cout << v2[j]<<" ";
std::cout << "\n";
std::cout << sign<<"\n";
std::cout <<find(electrons,v2, table)<<"\n";*/

         exit (0);
}


void av(int n, dcmplx *in, dcmplx *out,int tlen, Telt *T)
{
   int i, j;

   for (i=0; i<n; i++) out[i] = 0;
   for (i=0; i<tlen; i++) out[T[i].i] += in[T[i].j] * T[i].v;
}

int construct_hamiltonian( Site* lattice, int *table, unsigned long long
int *b,int electrons,int dimension, Telt *T)
{
     int i,state,j, k,l,m,n ,site,neighbor,ok,tlen,sign,*v1,*v2;
     double interaction;
     tlen=0;
     v1=(int*) malloc(electrons*sizeof(int));
     v2=(int*) malloc(electrons*sizeof(int));


     for (i=0;i<dimension;i++)
     {
     interaction=0;
     for (j=0;j<electrons;j++)
     {

         site=table[i*electrons+j]-1;
         for (k=0;k<COORD;k++)
         {
             neighbor=lattice[site].neighbors[k];

             ok=1;
             m=0;
             while (ok && (m<electrons))
             {
                 if ((table[i*electrons+m]-1)==neighbor) ok=0;
                 m++;
             }

             if (ok)
             {
                 memcpy(v1,&table[i*electrons], electrons*sizeof(int));



                 sign=hop(electrons, v1, v2,j, neighbor+1);
                 if (sign==0) sign=1; else sign=-1;

                 state=find(dimension,electrons,v2, table, b);

                 T[tlen].i=i;
                 T[tlen].j=state;
                 T[tlen].v=lattice[site].t[k];
                 T[tlen].v=T[tlen].v*(double)sign;
                 tlen++;


             }
             else
             {
             interaction=interaction+lattice[site].V[k];

             }
         }
     }
     if (interaction !=0)
     {
         T[tlen].i=i;
         T[tlen].j=i;
         T[tlen].v=interaction;
         tlen++;

     }
     }
     return tlen;

}


unsigned long long int C(int n, int k) //the number of combination of k
out of n
{

     int i;
     unsigned long long int c;

     c=1ULL;
     for (i=0;i<k;i++) c=c*(n-i);
     for (i=0;i<k;i++) c=c/(i+1);

     return c;

}

int combinations ( int n, int k, unsigned long long int *b,int *tab)
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


     for (i=1;i<C(n,k);i++)
     {
         //y= x | (x - 1);
         //x = (y + 1) | (((~y & -~y) - 1) >> (__builtin_ctz(x) + 1));


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



int find(int dimension,int k, int *v, int *T, unsigned long long int *b)
// find the position of a given combination v (vector of length k) in the table of all combinations T (array of length C(n,k)*k)
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
                {//cout<<"found at "<<mid<<"\n";
          return mid;}

     }
     //cout <<mid<<"\n";
}



int hop(int k, int *v1, int *v2,int n, int j) //given a state v1 (a combination corresponding to k occupied sites) generates the state v2
//obtained by hopping from site v1[n] to site j (which should not be in v1), and outputs the fermionic sign
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
