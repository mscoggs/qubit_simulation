#include <string.h>
#include <complex.h>

#include "linear_algebra.h"


void diag_hermitian_real_double(int N,  double *A, double *v_diag,double *e_vals){
	char JOBZ='V',UPLO='U';
	int LDA=N,LWORK=-1, INFO;
	double *WORK;
	WORK= new double[1]();

	memcpy(v_diag,A,N*N*sizeof(double));

	dsyev_(&JOBZ, &UPLO, &N, v_diag, &LDA, e_vals, WORK, &LWORK, &INFO );
	if (INFO !=0) printf("\n\n\nDIAGONALIZATION ERROR, INFO = %i\n\n\n", INFO);
 	LWORK=WORK[0];
 	delete[] WORK;
	WORK= new double[LWORK]();

 	dsyev_(&JOBZ, &UPLO, &N, v_diag, &LDA, e_vals, WORK, &LWORK, &INFO );
 	if (INFO !=0) printf("\n\n\nDIAGONALIZATION ERROR, INFO = %i\n\n\n", INFO);
	delete[] WORK;
}



void exp_diaganolized_real_matrix(double *hamiltonian, double *v_diag, double* e_vals, int N, double time_step){
	char TRANSA = 'N', TRANSB = 'N';
	int i, *IPIV, LWORK=N*N, INFO, LDA=N, LDB=N, LDC=N;
	double *exp_D, *temp_mat, *Vdag_z, *Vdag_z_inv, *WORK,ALPHA[2], BETA[2];
	ALPHA[0]=1.0, ALPHA[1]=0.0;
	BETA[0]=0.0, BETA[1]=0.0;

	IPIV = new int[N]();
	WORK =  new double[LWORK]();
	exp_D = new double[2*N*N]();
	temp_mat = new double[2*N*N]();
	Vdag_z = new double[2*N*N]();
	Vdag_z_inv = new double[2*N*N]();

	for (i =0; i<N*N*2; i++) hamiltonian[i]=0;
	for (i =0; i<N*N; i++) Vdag_z[2*i] = v_diag[i];

	for (i =0; i<N; i++) exp_D[2*(i*N+i)] = cos(-time_step*e_vals[i]), exp_D[2*(i*N+i)+1] = sin(-time_step*e_vals[i]);

	dgetrf_(&N,&N,v_diag,&N,IPIV,&INFO);
	dgetri_(&N,v_diag,&N,IPIV,WORK,&LWORK,&INFO);//inverting vdag

	for (i =0; i<N*N; i++) Vdag_z_inv[2*i] = v_diag[i];
	zgemm_(&TRANSA, &TRANSB, &N, &N, &N, ALPHA, Vdag_z, &LDA, exp_D, &LDB, BETA, temp_mat, &LDC); //matrix mult
	zgemm_(&TRANSA, &TRANSB, &N, &N, &N, ALPHA, temp_mat, &LDA, Vdag_z_inv, &LDB, BETA, hamiltonian, &LDC); //matrix mult

	delete[] exp_D, delete[] temp_mat, delete[] Vdag_z, delete[] Vdag_z_inv, delete[] IPIV, delete[] WORK;
}



void exp_complex_double_matrix_pade(int N, double *A, double *B){
	int M=N,K=N,ii,jj,kk,s,p,q,INFO,LDA=N,LDB=N,LDC=N,NRHS=N, *IPIV;
	double *row_norm,*X,*Y,*Z,*E,*D,norm,c, ALPHA[2], BETA[2];
	char TRANSA='N',TRANSB='N',TRANS='N';
	ALPHA[0]=1.0,ALPHA[1]=0.0;
	BETA[0]=0.0,BETA[1]=0.0;

	row_norm = new double[N]();
	X= new double[2*N*N]();
	Y= new double[2*N*N]();
	Z= new double[2*N*N]();
	E= new double[2*N*N]();
	D= new double[2*N*N]();
	IPIV= new int[N]();
	memcpy(Z,A,2*N*N*sizeof(double));

	for(ii=0;ii<N;ii++){
		row_norm[ii]=0.0;
		for (jj=0;jj<N;jj++) row_norm[ii]=row_norm[ii]+cabs(A[2*(jj+N*ii)]+I*A[2*(jj+N*ii)+1]);
	}

	norm=row_norm[0];
	for(ii=1;ii<N;ii++) if (row_norm[ii]>norm) norm=row_norm[ii];

	s=(int) floor(log2(norm)+2);
	if (s<0) s=0;

	for(ii=0;ii<2*N*N;ii++) Z[ii]=Z[ii]/pow(2.0,s);

	memcpy(X,Z,2*N*N*sizeof(double));
	c=0.5,p=1,q=6;

	for(ii=0;ii<2*N*N;ii++) E[ii]=c*Z[ii];
	for(ii=0;ii<2*N*N;ii++) D[ii]=-c*Z[ii];
	for(ii=0;ii<N;ii++) E[2*(ii+N*ii)]=E[2*(ii+N*ii)]+1.0;
	for(ii=0;ii<N;ii++) D[2*(ii+N*ii)]=D[2*(ii+N*ii)]+1.0;

	for (kk=2;kk<=q;kk++){
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
	if (INFO !=0) printf("\n\n\nERROR, INFO = %i\n\n\n", INFO);
	zgetrs_( &TRANS, &N, &NRHS,D,    &LDA, IPIV,E, &LDB, &INFO);//solving a system of equations
	if (INFO !=0) printf("\n\n\nERROR, INFO = %i\n\n\n", INFO);

	for (kk=1;kk<=s;kk++){
		memcpy(X,E,2*N*N*sizeof(double));
		zgemm_ (&TRANSA, &TRANSB, &M, &N, &K, ALPHA, E, &LDA, X, &LDB, BETA, Y, &LDC);//matrixmultiplication
		if (INFO !=0) printf("\n\n\nERROR, INFO = %i\n\n\n", INFO);
		memcpy(E,Y,2*N*N*sizeof(double));
	}
	memcpy(B,E,2*N*N*sizeof(double));

	delete[] row_norm, delete[] X, delete[] Y, delete[] Z, delete[] E, delete[] D, delete[] IPIV;
}



void matrix_vector_mult(double *matrix, double *state, int N){
	char TRANS = 'N';
	int i,INCX = 1,INCY = 1,LDA = N,M=N;
	double *result, ALPHA[2], BETA[2];
	ALPHA[0]=1.0,ALPHA[1]=0.0;
	BETA[0]=0.0,BETA[1]=0.0;

	result = new double[N*2]();

	zgemv_(&TRANS, &M, &N,ALPHA,matrix, &LDA, state, &INCX, BETA, result, &INCY);
	memcpy(state,result, 2*N*sizeof(double));

	delete[] result;
}
