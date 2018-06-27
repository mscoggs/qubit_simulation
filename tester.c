#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define COORD 14
#define CANT acos(1.0/3.0)/2.0
#define INTERACTION 0.35
#define nu 1/9
#define DISCRET 0
#define EIGENVALUES 40
#define T1 -0.13
#define T2 -0.02
#define T3 0.05







unsigned long long int choose(int n, int k) //the number of combination of k out of n
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

     for (i=1;i<choose(n,k);i++)
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






int construct_hamiltonian( Site* lattice, int *table, unsigned long long int *b,int electrons,int dimension, Telt *T)
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


int main (int argc, char *argv[])
{
     int *ising,*table,*v, *v2,i,j,k,sign, nev,tlen,Nx,Ny,sites, electrons, dimension;
     unsigned long long int *b;
     FILE *fp;
     Site *lattice;
     Telt *T;
     dcmplx **Evecs, *Evals;

     //Nx = atoi(argv[1]);
     //Ny = atoi(argv[2]);

     sites=9;
     electrons=sites*nu; // 1
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

     tlen=construct_hamiltonian( lattice,table,b, electrons, dimension, T);
     exit (0);
}

'''
int mainORIGINALORIGINALOGOGOGOGOGOOGOG(int argc, char *argv[])
{
     int *ising,*table,*v, *v2,i,j,k,sign, nev,tlen,Nx,Ny,sites, electrons, dimension;
     unsigned long long int *b;
     FILE *fp;
     Site *lattice;
     Telt *T;
     dcmplx **Evecs, *Evals;
     double phi;

     Nx = atoi(argv[1]);
     Ny = atoi(argv[2]);
     phi= atof(argv[3]);

     sites=3*Nx*Ny;

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

     tlen=construct_hamiltonian( lattice,table,b, electrons, dimension, T);
     exit (0);
}''
