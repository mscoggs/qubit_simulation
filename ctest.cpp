#include <iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <cblas.h>


using namespace std;
#define nx 4

//https://en.wikipedia.org/wiki/Pad%C3%A9_approximant
//http://www.maths.manchester.ac.uk/~higham/talks/exp09.pdf


void diag_hermitian_real_double(int N,  double *A, double *Vdag,double *D)
{


     char JOBZ,UPLO;
     int LDA,LWORK;
     double *WORK;
     int INFO;
     JOBZ='V';
     UPLO='U';
     LDA=N;
     LWORK=-1;
     WORK=malloc(sizeof(double));
     memcpy(Vdag,A,N*N*sizeof(double));
     dsyev_(&JOBZ, &UPLO, &N, Vdag, &LDA, D, WORK, &LWORK, &INFO );
     LWORK=WORK[0];
     free(WORK);
     WORK=malloc(LWORK*sizeof(double));
     dsyev_(&JOBZ, &UPLO, &N, Vdag, &LDA, D, WORK, &LWORK, &INFO );
     if (INFO !=0) printf("DIAGONALIZATION ERROE\n");
     free(WORK);
}

  int main()
  {

  }
  free(IPIV);
}
