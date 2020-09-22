#ifndef __LINALG_H_INCLUDED__
#define __LINALG_H_INCLUDED__


extern "C" int zgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K, double *ALPHA, double *Z, int *LDA, double *X, int *LDB, double *BETA, double *Y, int *LDC); //complex matrix*matrix mult, odd indices hold imaginary values.
extern "C" int zgemv_(char *TRANS, int *M, int *N,double *ALPHA,double *A, int *LDA, double *X, int *INCX, double *BETA, double *Y, int *INCY); //complex matrix-vector mult, odd indices hold imaginary values.
extern "C" int dsyev_(char *JOBZ, char *UPLO, int *N, double *v_diag, int *LDA, double *e_vals, double *WORK, int *LWORK, int *INFO);//diagonalization, returns the eigenvectors in v_diag and eigenvalues in e_vals.
extern "C" int zgetrf_ (int *M, int *N, double *D, int *LDA, int *IPIV, int *INFO);//A matrix factorization
extern "C" int zgetrs_( char *TRANS, int *N, int *NRHS,double *D, int *LDA, int *IPIV,double *E, int *LDB, int *INFO);//solves a system of equations
extern "C" void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
extern "C" void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
extern "C" int dgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K, double *ALPHA, double *Z, int *LDA, double *X, int *LDB, double *BETA, double *Y, int *LDC); //real values matrix*matrix mult
extern "C" double zdotc_(int *N, double*ZX,int *INCX, double *ZY, int *INCY);//dots the complex conjugate of ZX with ZY


/**
    Diagonalizing a real square matrix A. Stores the eigenvectors in v_diag and the eigenvalues in e_vals

    @param N the dimension of the matrix
    @param A the matrix to be diagonalized
    @param v_diag an array (size N*N) that will contain the eigen vectors of the matrix A
    @param e_vals an array (size N) that will be filled with the eigen values of the matrix A
*/
void diag_hermitian_real_double(int N,  double *A, double *v_diag,double *e_vals);


/**
    Calculating the exponential of a diagonalized decomposition where the matrix A = v_diag*e_vals*Vdag_inv.
        calculating v_diag*exp(e_vals)*Vdag_inv =exp(A), storing the result in hamiltonian

    @param ham_real the real matrix that will hold the result
    @param v_diag an array (size N*N) that contains the eigen vectors of the matrix A
    @param e_vals an array (size N) that contains the eigen values of the matrix A
    @param N the dimension of the matrix
    @param time_step the scalar that multiplies the matrix before being exponentiated
*/
void exp_diaganolized_real_matrix(double *ham_real, double *v_diag, double* e_vals, int N, double time_step);


/**
    Using the Pade Approximation to calulate exp(A)

    @param N the dimension of the matrix
    @param A a complex matrix storing real values at
        even indices and their imaginary parts at the odd indices
    @param B this holds the result
*/
void exp_complex_double_matrix_pade(int N, double *A, double *B);


/**
    Matrix vector multiplication, state=matrix*state

    @param matrix the N*N matrix
    @param state  the vector of size N
    @param N the dimension of the matrix
*/
void matrix_vector_mult(double *matrix, double *state, int N);

#endif
