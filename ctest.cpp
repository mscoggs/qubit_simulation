#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <complex>
#include <math.h>

void matrix_vector_mult_c(double *exp_matrix, double *psi, double *psi_result,int dimension);


using namespace std;
#define Ncount 20
void printArray(int scalar, int offset);

double list2[][2]= {{1,2},{1,2}};

double JList[18]=
{7.0,8.1,9.2,10.3,11.4,12.5,13.6,14.7,15.8,16.9,18.0,19.1,20.2,21.3,22.4,23.5,24.6,25.7};
double newList[2];
  int main()
  {
    double *newList;
    newList = (double*) malloc(sizeof(double)*2);
    int row = 1;
    int col =1;
    newList[0]=row;
    newList[1]=col;
    printf("%d\n",newList[1]);
    printf("%i\n",(uintptr_t)JList[4]);
  }

  void printArray(int scalar, int offset){
    int i,j;
    printf("{");
    for (i=0;i<Ncount;i++)
    {
      printf("{%.1f",i*1.1+offset);
      for (j=1;j<3*scalar;j++) printf(",%.1f",((i+j)*1.1+offset));
      if (i<Ncount-1)
      {
        if ((i%4)==3) printf("},\n");
        else printf("},");
      }
      else printf("}};\n\n");
    }
  }






  void matrix_vector_mult_c(double *exp_matrix, double *psi,double *psi_result, int dimension)
  {
    int i,j;
    double *psi_temp;
    psi_temp = (double*) malloc (2*dimension*sizeof(double));
    for(i=0; i<dimension; i++)
    {
      double resultR = 0.0;
      double resultI = 0.0;

      for (j=0; j<dimension; j++)
      {
        resultR += exp_matrix[2*j+i*dimension]*psi[2*j];
        resultI += exp_matrix[2*j+i*dimension+1]*psi[2*j+1];
      }
      psi_temp[2*i] = resultR;
      psi_temp[2*i+1] = resultI;
    }
    memcpy(psi_result, psi_temp, 2*dimension*sizeof(double));
    for (i=0;i<dimension*2;i++) printf("%f\n", psi_result[i]);

  }
