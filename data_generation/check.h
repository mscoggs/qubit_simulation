#ifndef __CHECK_H_INCLUDED__
#define __CHECK_H_INCLUDED__


/**
    Checks the commutator of two complex hamiltonians. This reports an error if the commutator is zero, which would result
        in no change during evolution

    @param N the dimension of the quantum system and the hamiltonians
    @param A one the the hamiltonians of size 2*N X 2*N, where odd indices hold imaginary components
    @param B one the the hamiltonians of size 2*N X 2*N

    @return true if they do commute, false if they don't
*/
bool check_commutator(int N, double* A, double* B);


/**
    Confirms that the hamiltonian is unitary, if not it reports an error

    @param hamiltonian the hamiltonian that will be checked, of size (2*N by 2*N)
    @param N the dimension of the quantum system and the hamiltonians
*/
void check_unitary(double* hamiltonian, int N);


/**
    Confirms that the hamiltonian is hermitian, if not it reports an error

    @param hamiltonian the hamiltonian that will be checked, of size (2*N by 2*N)
    @param N the dimension of the quantum system and the hamiltonians
*/
void check_hermicity(double* hamiltonian, int N);


/**
    Confirms that the square of the weights for a given state add up to 1. For these quantum state, state = sum {|c_i ^2|state_i}
        where state_i is the ith eigenvector of the hamiltonian and sum{c_i ^ 2} should be 1.

    @param state the state we're using the check the weights
    @param hamiltonian the hamiltonian that will be checked, of size (2*N by 2*N)
    @param N the dimension of the quantum system and the hamiltonians
*/
void check_weights(double* state, double* hamiltonian, int N);


/**
    Confirms that the norm of the state is 1.

    @param state the state we're checking, of size (2*N)
    @param N the dimension of the quantum system and the hamiltonians
*/
void check_norm(double* state, int N);

#endif
