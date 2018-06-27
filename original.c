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






void construct_lattice( Site* lattice, int Nx,int Ny, double V,double
delta, int *ising, double phi, double tp, double tpp, double tpph)
{
int sites,ii,jj,kk,aa,bb;
sites=Nx*Ny*3;
double beta[sites], alpha[sites];



for (ii=0;ii<Nx;ii++)for (jj=0;jj<Ny;jj++)
{
alpha[3*(Nx*jj+ii)+0]=M_PI/6;
alpha[3*(Nx*jj+ii)+1]=5*M_PI/6;
alpha[3*(Nx*jj+ii)+2]=M_PI/2;
beta[3*(Nx*jj+ii)+0]=M_PI/2-delta;
beta[3*(Nx*jj+ii)+1]=M_PI/2-delta;
beta[3*(Nx*jj+ii)+2]=M_PI/2+delta;
}
for (ii=0;ii<Nx;ii++)for (jj=0;jj<Ny;jj++)for (kk=0;kk<3;kk++)
if ((ising[3*(Nx*jj+ii)+kk])<0)
{
alpha[3*(Nx*jj+ii)+kk]=M_PI+alpha[3*(Nx*jj+ii)+kk];
beta[3*(Nx*jj+ii)+kk]=M_PI-beta[3*(Nx*jj+ii)+kk];
}






//rows hopping
for (jj=0;jj<Ny;jj++)for (ii=0;ii<Nx;ii++)
{
aa=3*(Nx*jj+ii)+1;
bb=3*(Nx*jj+(ii+1)%Nx);
lattice[aa].neighbors[3]=bb;
lattice[bb].neighbors[2]=aa;
lattice[aa].t[3]=(cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))+I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb]));
lattice[bb].t[2]=(cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))-I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb]));
lattice[aa].t[3]=lattice[aa].t[3]*(cos(phi/2/Nx)+I*sin(phi/2/Nx));
lattice[bb].t[2]=lattice[bb].t[2]*(cos(phi/2/Nx)-I*sin(phi/2/Nx));

lattice[aa].V[3]=V/2;
lattice[bb].V[2]=V/2;

/*A[2*(bb+sites*aa)]=(cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]));
A[2*(aa+sites*bb)]=(cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]));
A[2*(bb+sites*aa)+1]=(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb]));
A[2*(aa+sites*bb)+1]=-(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb]));*/
}
//column hopping
for (jj=0;jj<Ny;jj++)for (ii=0;ii<Nx;ii++)
{
aa=3*(Nx*jj+ii)+2;
bb=3*(Nx*((jj+1)%Ny)+ii);
lattice[aa].neighbors[2]=bb;
lattice[bb].neighbors[3]=aa;
lattice[aa].t[2]=(cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))+I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb]));
lattice[bb].t[3]=(cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))-I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb]));

lattice[aa].V[2]=V/2;
lattice[bb].V[3]=V/2;



/*A[2*(bb+sites*aa)]=(cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]));
A[2*(aa+sites*bb)]=(cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]));
A[2*(bb+sites*aa)+1]=(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb]));
A[2*(aa+sites*bb)+1]=-(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb]));*/
}
//diagonal hopping
for (jj=0;jj<Ny;jj++)for (ii=0;ii<Nx;ii++)
{
aa=3*(Nx*jj+ii)+2;
bb=3*(Nx*((jj+1)%Ny)+(Nx+ii-1)%Nx)+1;
lattice[aa].neighbors[3]=bb;
lattice[bb].neighbors[2]=aa;
lattice[aa].t[3]=(cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))+I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb]));
lattice[bb].t[2]=(cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))-I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb]));
lattice[aa].t[3]=lattice[aa].t[3]*(cos(phi/2/Nx)-I*sin(phi/2/Nx));
lattice[bb].t[2]=lattice[bb].t[2]*(cos(phi/2/Nx)+I*sin(phi/2/Nx));

lattice[aa].V[3]=V/2;
lattice[bb].V[2]=V/2;

}
//internal hopping
for (jj=0;jj<Ny;jj++)for (ii=0;ii<Nx;ii++)
{
aa=3*(Nx*jj+ii);
bb=3*(Nx*jj+ii)+1;
lattice[aa].neighbors[0]=bb;
lattice[bb].neighbors[1]=aa;
lattice[aa].t[0]=(cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))+I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb]));
lattice[bb].t[1]=(cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))-I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb]));
lattice[aa].t[0]=lattice[aa].t[0]*(cos(phi/2/Nx)+I*sin(phi/2/Nx));
lattice[bb].t[1]=lattice[bb].t[1]*(cos(phi/2/Nx)-I*sin(phi/2/Nx));


lattice[aa].V[0]=V/2;
lattice[bb].V[1]=V/2;


aa=3*(Nx*jj+ii);
bb=3*(Nx*jj+ii)+2;
lattice[aa].neighbors[1]=bb;
lattice[bb].neighbors[0]=aa;
lattice[aa].t[1]=(cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))+I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb]));
lattice[bb].t[0]=(cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))-I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb]));

lattice[aa].V[1]=V/2;
lattice[bb].V[0]=V/2;

aa=3*(Nx*jj+ii)+1;
bb=3*(Nx*jj+ii)+2;
lattice[aa].neighbors[0]=bb;
lattice[bb].neighbors[1]=aa;
lattice[aa].t[0]=(cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))+I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb]));
lattice[bb].t[1]=(cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))-I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb]));
lattice[aa].t[0]=lattice[aa].t[0]*(cos(phi/2/Nx)-I*sin(phi/2/Nx));
lattice[bb].t[1]=lattice[bb].t[1]*(cos(phi/2/Nx)+I*sin(phi/2/Nx));

lattice[aa].V[0]=V/2;
lattice[bb].V[1]=V/2;
}

//NNN hopping
for (jj=0;jj<Ny;jj++)for (ii=0;ii<Nx;ii++)
{


aa=3*(Nx*jj+ii);
bb=3*(Nx*((jj+1)%Ny)+(ii-1+Nx)%Nx)+1;


lattice[aa].neighbors[4]=bb;
lattice[bb].neighbors[5]=aa;
lattice[aa].t[4]=((cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))+I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb])))*tp;
lattice[bb].t[5]=((cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))-I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb])))*tp;
lattice[aa].t[4]=lattice[aa].t[4]*(cos(phi/Nx/2)-I*sin(phi/Nx/2));
lattice[bb].t[5]=lattice[bb].t[5]*(cos(phi/Nx/2)+I*sin(phi/Nx/2));

lattice[aa].V[4]=0.0;
lattice[bb].V[5]=0.0;



aa=3*(Nx*jj+ii);
bb=3*(Nx*((jj+1)%Ny)+(ii-1+Nx)%Nx);


lattice[aa].neighbors[8]=bb;
lattice[bb].neighbors[9]=aa;
lattice[aa].t[8]=((cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))+I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb])))*tpph;
lattice[bb].t[9]=((cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))-I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb])))*tpph;
lattice[aa].t[8]=lattice[aa].t[8]*(cos(phi/Nx)-I*sin(phi/Nx));
lattice[bb].t[9]=lattice[bb].t[9]*(cos(phi/Nx)+I*sin(phi/Nx));

lattice[aa].V[8]=0.0;
lattice[bb].V[9]=0.0;



aa=3*(Nx*jj+ii);
bb=3*(Nx*jj+(ii-1+Nx)%Nx)+2;


lattice[aa].neighbors[5]=bb;
lattice[bb].neighbors[7]=aa;
lattice[aa].t[5]=((cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))+I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb])))*tp;
lattice[bb].t[7]=((cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))-I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb])))*tp;
lattice[aa].t[5]=lattice[aa].t[5]*(cos(phi/Nx)-I*sin(phi/Nx));
lattice[bb].t[7]=lattice[bb].t[7]*(cos(phi/Nx)+I*sin(phi/Nx));

lattice[aa].V[5]=0.0;
lattice[bb].V[4]=0.0;


aa=3*(Nx*jj+ii);
bb=3*(Nx*((jj+Ny-1)%Ny)+(ii+1)%Nx)+2;


lattice[aa].neighbors[6]=bb;
lattice[bb].neighbors[4]=aa;
lattice[aa].t[6]=((cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))+I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb])))*tp;
lattice[bb].t[4]=((cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))-I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb])))*tp;
lattice[aa].t[6]=lattice[aa].t[6]*(cos(phi/Nx)+I*sin(phi/Nx));
lattice[bb].t[4]=lattice[bb].t[4]*(cos(phi/Nx)-I*sin(phi/Nx));

lattice[aa].V[6]=0.0;
lattice[bb].V[4]=0.0;




aa=3*(Nx*jj+ii);
bb=3*(Nx*((jj+Ny-1)%Ny)+ii)+1;


lattice[aa].neighbors[7]=bb;
lattice[bb].neighbors[6]=aa;
lattice[aa].t[7]=((cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))+I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb])))*tp;
lattice[bb].t[6]=((cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))-I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb])))*tp;
lattice[aa].t[7]=lattice[aa].t[7]*(cos(phi/Nx/2)+I*sin(phi/Nx/2));
lattice[bb].t[6]=lattice[bb].t[6]*(cos(phi/Nx/2)-I*sin(phi/Nx/2));

lattice[aa].V[7]=0.0;
lattice[bb].V[6]=0.0;


aa=3*(Nx*jj+ii);
bb=3*(Nx*jj+(ii+1)%Nx);


lattice[aa].neighbors[10]=bb;
lattice[bb].neighbors[12]=aa;
lattice[aa].t[10]=((cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))+I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb])))*tpp;
lattice[bb].t[12]=((cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))-I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb])))*tpp;
lattice[aa].t[10]=lattice[aa].t[10]*(cos(phi/Nx)+I*sin(phi/Nx));
lattice[bb].t[12]=lattice[bb].t[12]*(cos(phi/Nx)-I*sin(phi/Nx));

lattice[aa].V[10]=0.0;
lattice[bb].V[12]=0.0;

aa=3*(Nx*jj+ii);
bb=3*(Nx*((jj+1)%Ny)+ii);


lattice[aa].neighbors[11]=bb;
lattice[bb].neighbors[13]=aa;
lattice[aa].t[11]=((cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))+I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb])))*tpp;
lattice[bb].t[13]=((cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))-I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb])))*tpp;


lattice[aa].V[11]=0.0;
lattice[bb].V[13]=0.0;


aa=3*(Nx*jj+ii)+1;
bb=3*(Nx*jj+(ii+1)%Nx)+1;


lattice[aa].neighbors[13]=bb;
lattice[bb].neighbors[11]=aa;
lattice[aa].t[13]=((cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))+I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb])))*tpp;
lattice[bb].t[11]=((cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))-I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb])))*tpp;
lattice[aa].t[13]=lattice[aa].t[13]*(cos(phi/Nx)+I*sin(phi/Nx));
lattice[bb].t[11]=lattice[bb].t[11]*(cos(phi/Nx)-I*sin(phi/Nx));

lattice[aa].V[13]=0.0;
lattice[bb].V[11]=0.0;




aa=3*(Nx*jj+ii)+1;
bb=3*( Nx*( (jj+1)%Ny )+ (ii-1+Nx)%Nx )+1;


lattice[aa].neighbors[10]=bb;
lattice[bb].neighbors[12]=aa;
lattice[aa].t[10]=((cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))+I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb])))*tpp;
lattice[bb].t[12]=((cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))-I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb])))*tpp;
lattice[aa].t[10]=lattice[aa].t[10]*(cos(phi/Nx)-I*sin(phi/Nx));
lattice[bb].t[12]=lattice[bb].t[12]*(cos(phi/Nx)+I*sin(phi/Nx));

lattice[aa].V[10]=0.0;
lattice[bb].V[12]=0.0;

aa=3*(Nx*jj+ii)+2;
bb=3*(Nx*((jj+1)%Ny)+ii)+2;


lattice[aa].neighbors[10]=bb;
lattice[bb].neighbors[12]=aa;
lattice[aa].t[10]=((cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))+I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb])))*tpp;
lattice[bb].t[12]=((cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))-I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb])))*tpp;


lattice[aa].V[10]=0.0;
lattice[bb].V[12]=0.0;



aa=3*(Nx*jj+ii)+2;
bb=3*(Nx*((jj+1)%Ny)+(ii-1+Nx)%Nx)+2;


lattice[aa].neighbors[11]=bb;
lattice[bb].neighbors[13]=aa;
lattice[aa].t[11]=((cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))+I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb])))*tpp;
lattice[bb].t[13]=((cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))-I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb])))*tpp;
lattice[aa].t[11]=lattice[aa].t[11]*(cos(phi/Nx)-I*sin(phi/Nx));
lattice[bb].t[13]=lattice[bb].t[13]*(cos(phi/Nx)+I*sin(phi/Nx));

lattice[aa].V[11]=0.0;
lattice[bb].V[13]=0.0;




aa=3*(Nx*jj+ii)+1;
bb=3*(Nx*((jj+1)%Ny)+ii)+1;


lattice[aa].neighbors[9]=bb;
lattice[bb].neighbors[8]=aa;
lattice[aa].t[9]=((cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))+I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb])))*tpph;
lattice[bb].t[8]=((cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))-I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb])))*tpph;


lattice[aa].V[9]=0.0;
lattice[bb].V[8]=0.0;


aa=3*(Nx*jj+ii)+1;
bb=3*(Nx*jj+(ii+1)%Nx)+2;


lattice[aa].neighbors[7]=bb;
lattice[bb].neighbors[5]=aa;
lattice[aa].t[7]=((cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))+I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb])))*tp;
lattice[bb].t[5]=((cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))-I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb])))*tp;
lattice[aa].t[7]=lattice[aa].t[7]*(cos(phi/Nx/2)+I*sin(phi/Nx/2));
lattice[bb].t[5]=lattice[bb].t[5]*(cos(phi/Nx/2)-I*sin(phi/Nx/2));

lattice[aa].V[7]=0.0;
lattice[bb].V[5]=0.0;



aa=3*(Nx*jj+ii)+1;
bb=3*(Nx*((jj+Ny-1)%Ny)+ii)+2;


lattice[aa].neighbors[4]=bb;
lattice[bb].neighbors[6]=aa;
lattice[aa].t[4]=((cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))+I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb])))*tp;
lattice[bb].t[6]=((cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))-I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb])))*tp;
lattice[aa].t[4]=lattice[aa].t[4]*(cos(phi/Nx/2)-I*sin(phi/Nx/2));
lattice[bb].t[6]=lattice[bb].t[6]*(cos(phi/Nx/2)+I*sin(phi/Nx/2));

lattice[aa].V[4]=0.0;
lattice[bb].V[6]=0.0;





aa=3*(Nx*jj+ii)+2;
bb=3*(Nx*jj+(ii+1)%Nx)+2;


lattice[aa].neighbors[9]=bb;
lattice[bb].neighbors[8]=aa;
lattice[aa].t[9]=((cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))+I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb])))*tpph;
lattice[bb].t[8]=((cos(beta[aa]/2)*cos(beta[bb]/2)+sin(beta[aa]/2)*sin(beta[bb]/2)*cos(alpha[aa]-alpha[bb]))-I*(sin(beta[aa]/2)*sin(beta[bb]/2)*sin(alpha[aa]-alpha[bb])))*tpph;
lattice[aa].t[9]=lattice[aa].t[9]*(cos(phi/Nx)+I*sin(phi/Nx));
lattice[bb].t[8]=lattice[bb].t[8]*(cos(phi/Nx)-I*sin(phi/Nx));

lattice[aa].V[9]=0.0;
lattice[bb].V[8]=0.0;


}



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
