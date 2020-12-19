#ifndef __OPERATIONS_H_INCLUDED__
#define __OPERATIONS_H_INCLUDED__

#include <gsl/gsl_rng.h>

#include "parameters.h"


/**
		Return the ground state of a hamiltonian

		@param N the dimension of our system
		@param hamiltonian the complex hamiltonian of dimension 2*Nx2*N
		@param ground_state will hold the desired ground state, of size 2*N
*/
double get_ground_E(int N, double *hamiltonian);


/**
		Returns the ground state energy of a hamiltonian

		@param N the dimension of our system
		@param hamiltonian the complex hamiltonian of dimension 2*Nx2*N

		@return min_E the ground state energy
*/
void get_ground_state(int N, double *hamiltonian, double* ground_state);


/**
		Return the expectation value between an evolved state state and the target hamiltonian , expressed as <state|hamiltonian|state> if USE_ENERGY_DISTANCE is true,
    otherwise it returns the complex_dot_squared between the target state and the best state after evolution (state)

		@param N the dimension of our system
		@param state the evolved state
		@param hamiltonian the complex hamiltonian of dimension 2*Nx2*N

		@return result the expectation value of <state|H|state> OR the complex_dot_squared
*/
double cost(double *target_state, int N, double *state, double *ham_target, bool energy = USE_ENERGY_DISTANCE);


/**
		Given a state b generates the state v obtained by hopping from site n to site m (which should not be in b), and outputs the fermionic sign

		@param b the start state
		@param v the target state
		@param n the occupied site we start from
		@param m the unoccupied site we hop to

		@return z_count%2 fermionic sign
*/
int hop(unsigned long long int b, unsigned long long int *v,int n, int j);


/**
		Find the position of a given state v in the collection of all states b

		@param N the dimension of the system
		@param v the state whos index we're looking for
		@param b a collection of all binary-represented states

		@return mid+1 the desired index
*/
int find(int N,unsigned long long int *v, unsigned long long int *b);


/**
		A simple choose algorithm, NUMBER_OF_SITES choose num_occupants, giving the dimension of our system

		@param num_occupants the number to be chosen, the number of occupied sites on our lattice

		@return c the NUMBER_OF_SITES choose num_occupants, or (NUMBER_OF_SITES! / num_occupants!)
*/
unsigned long long int choose(int num_occupants);


/**
		Getting the binary representation of all our states (where a 1 means occupied, and 0 means unoccupied) and storing them in b. Collecting the corresponding site numbers for a site and storing it in tab.
		Ex: for a 2x2 lattice, with 2 occupants, N = 6 (4choose2). We'd have 6 different states, 2 occupants in those 4 locations. b[0]= 3 = 0011, meaning site 1 and 2 are occupied with 3 and 4 empty. This would
		correspond to tab[0] = 1, tab[1] = 2, meaning site 1 and 2 are occupied in the first state.

		@param N the dimension of our system
		@param num_occupants
		@param b the array holding all our states in binary representation
		@param tab the array holding the site which are occupied
*/
int combinations( int N, int num_occupants, unsigned long long int *b,int *tab);


/**
		Randomly generating a number in [lower, upper), including lower, up to upper (but not including)

		@param lower lower bound
		@param upper upper bound
		@param rng the rng key

		@return lower_to_upper the random number within our bounds
*/
double get_random_double(double lower, double upper, gsl_rng *rng);


/**
	  Creating a normalized distance from the initial energy to the target state energy

		@param initial the expectation value between the initial state and the target hamiltonian
		@param target the ground state of the target hamiltonian
		@param current the best energy post monte-carlo simulation

		@return distance the distance, in [0,1]
*/
double calc_distance(Simulation_Parameters& sim_params);


/**
    Calculates the tau for the next simulation

    @param sim_params contains all of the variables for the simulation
*/
void calc_tau(Simulation_Parameters& sim_params);

double calc_new_temperature(Simulation_Parameters& sim_params, int temp_index);

/**
		Build a NX*NY lattice, with sites 1 through NX*NY listed in a snaking pattern

		@param lattice the empty lattice to be filled with site numbers
*/
void construct_lattice(int lattice[NX][NY]);


/**
		Attaching a bond number for each neighbor pair, allowing assignment of a specific bond-value between two neighbors, where this value is stored in a NUMBER_OF_SITES by NUMBER_OF_SITES list

		@param bonds the NUMBER_OF_SITES by NUMBER_OF_SITES array that will hold the bond numbers
		@param lattice the lattice that holds the site numbers
*/
void assign_bonds(int *bonds, int lattice[NX][NY]);


/**
		Gets the neighbors of a site, returning all 4 neighbors if open_boundry coditions are true, otherwise just closed-boudnary neighbors

		@param site the site which we'll get the neighbors for
		@param neighbors the array which will hold the site numbers of the corresponding neighbors
		@param lattice the lattice that holds the site numbers

		@return count the total number of neighbors
*/
int get_neighbors(int site, int *neighbors, int lattice[NX][NY]);


/**
	  Updates the distances after a number of (NUM_SEEDS) monte carlo simulations. Checks to see if there's a convergence issue, meaning the most recent simulation was further away
    from the ground state than the previous. If this is the case, return false. If not, return true with the updates distances.

		@param sim_params an object that contains all of the parameters for the simulation

		@return bool true if the most recent simulation made progress towards the ground state, false if not
*/
bool update_distances(Simulation_Parameters& sim_params);


/**
    for all seeds, get the best energy, best evolved state, and best overlap squared between the evolved state and target state

    @param sim_params an object that contains all of the parameters for the simulation
*/
void get_best_seed(Simulation_Parameters& sim_params);


/**
    calculates the inner product between two states 

    @param size total number of complex numbers in each state, half the indices
    @param v1 the first state
    @param v2 the second state

    return total
*/
double complex_dot_squared(int size, double *v1, double *v2);


/**
    normlaizing the state by calculating its inner product and dividing each element by that
*/
void normalize_state(double *state, int N);


/**
    checking a few conditions, deciding if we should stop the MC simulation 
*/
bool exit_simulation(Simulation_Parameters& sim_params);

#endif
