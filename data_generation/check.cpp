#include <string.h>
#include <math.h>
#include <string>

#include "check.h"
#include "linear_algebra.h"
#include "print.h"
#include "parameters.h"



void check_commutator(int N, double* A, double* B){
	char TRANSA = 'N', TRANSB = 'N';
	int i, *IPIV, LWORK=N*N, INFO, LDA=N, LDB=N, LDC=N;
	double *C, *AB, *BA ,ALPHA[2], BETA[2];
	bool commute = true;
	ALPHA[0]=1.0, ALPHA[1]=0.0;
	BETA[0]=0.0, BETA[1]=0.0;


	C = new double[2*N*N]();
	BA = new double[2*N*N]();
	AB = new double[2*N*N]();

	zgemm_(&TRANSA, &TRANSB, &N, &N, &N, ALPHA, A, &LDA, B, &LDB, BETA, AB, &LDC); //matrix mult
	zgemm_(&TRANSA, &TRANSB, &N, &N, &N, ALPHA, B, &LDA, A, &LDB, BETA, BA, &LDC); //matrix mult
	for (i =0; i<N*N*2; i++) C[i] = AB[i] - BA[i];
	for (i =0; i<N*N*2; i++) if(C[i] < -0.001 || 0.001 < C[i]) commute = false;
	if(commute) printf("\n\n\nWARNING: THE TARGET AND INITIAL HAMILTONIAN COMMUTE\n\n\n");
	if(PRINT_COMMUTATOR) print_hamiltonian_complex(C, N);

	delete[] C, delete[] AB, delete[] BA;
}



void check_unitary(double* hamiltonian, int N){
	int i,*IPIV, LWORK=N*N, INFO, LDA=N, LDB=N, LDC=N;
	double *ham_t,*unitary, *WORK,ALPHA[2], BETA[2];
	char TRANSA = 'C', TRANSB = 'N';


	ALPHA[0]=1.0, ALPHA[1]=0.0;
	BETA[0]=0.0, BETA[1]=0.0;

	ham_t = new double[2*N*N]();
	unitary = new double[2*N*N]();
	memcpy(ham_t, hamiltonian, sizeof(double)*2*N*N);
	zgemm_(&TRANSA, &TRANSB, &N, &N, &N, ALPHA, ham_t, &LDA, hamiltonian, &LDB, BETA, unitary, &LDC); //matrix mult

	for(i=0;i<N;i++) unitary[2*(i*N+i)] = unitary[2*(i*N+i)] -1;
	for(i=0;i<N*N*2;i++) if(unitary[i] < -0.00000001 or unitary[i] > 0.00000001) printf("\n\n\nERROR, NON UNITARY ELEMENTS AT %i, VALUE: %f\n\n\n", i, unitary[i]), exit(0);
	delete[] ham_t, delete[] unitary;
}



void check_hermicity(double* hamiltonian, int N){
	int i,j;
	for(i=0;i<N;i++) for(j=0;j<N;j++) if(abs(hamiltonian[2*(j*N+i)] - hamiltonian[2*(i*N+j)]) > 0.00001) printf("\n\n\nERROR, NON HERMITIAN ELEMENT AT i,j: %i,%i\nElement_ij = %10.7f\nElement_ji = %10.7f\n\n\n", i*N+j, j*N+i, hamiltonian[2*(i*N+j)], hamiltonian[2*(j*N+i)]);
}



void check_weights(double* state, double* hamiltonian, int N){
	int i,j;
	double *evals, *v_diag, *ham_real, sum_real, sum_im, c_squared_sum=0;
	v_diag = new double[N*N]();
	evals = new double[N]();
	ham_real = new double[N*N]();

	for (i=0;i<N*N;i++) ham_real[i] = hamiltonian[2*i];

	diag_hermitian_real_double(N, ham_real,v_diag, evals);

	for(j=0;j<N;j++){
		sum_real=0, sum_im=0;
		for(i=0; i<N; i++) sum_real += state[2*i]*v_diag[j*N+i], sum_im += -state[2*i+1]*v_diag[j*N+i];
		c_squared_sum += sum_real*sum_real+sum_im*sum_im;
	}

	if(c_squared_sum>1.00001 or c_squared_sum <0.999999) printf("\n\n\nERROR, BAD WEIGHTS\n\n\n");

	delete[] v_diag, delete[] evals, delete[] ham_real;

}



void check_norm(double* state, int N){

	int i;
	double sum=0;
	for(i=0;i<N*2; i+=2) sum+= state[i]*state[i]+(state[i+1]*state[i+1]);
	if(sum>1.000000001 or sum<0.99999999) printf("\n\n\nNORM ERROR, SIZE: %f\n\n\n", sum),print_state(state, N), exit(0);
}
