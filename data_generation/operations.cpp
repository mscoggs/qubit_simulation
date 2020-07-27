#include <string.h>

#include "parameters.h"
#include "print.h"
#include "check.h"
#include "linear_algebra.h"
#include "operations.h"



void get_ground_state(int N, double *hamiltonian, double *ground_state){
	int i, min_i;
	double *evals, *v_diag, *ham_real, min_E=1000;

	evals = new double[N]();
	v_diag = new double[N*N]();
	ham_real = new double[N*N]();

	for(i=0;i<N*N;i++) ham_real[i] = hamiltonian[2*i];
	for(i=0; i<2*N; i++) ground_state[i] = 0.0;

	diag_hermitian_real_double(N, ham_real,v_diag, evals);

	for(i=0; i<N; i++) if(evals[i]<min_E) min_E = evals[i], min_i = i;
	for(i=0; i<N; i++) ground_state[i*2] = v_diag[(min_i)*N+i];

	delete[] evals, delete[] v_diag, delete[] ham_real;
}



double get_ground_E(int N, double *hamiltonian){
	int i;
	double *evals, *v_diag, *ham_real, min_E=1000;
	evals = new double[N]();
	v_diag = new double[N*N]();
	ham_real = new double[N*N]();

	for(i=0;i<N*N;i++) ham_real[i] = hamiltonian[2*i];

	diag_hermitian_real_double(N, ham_real,v_diag, evals);

	for(i=0; i<N; i++) if(evals[i]<min_E) min_E = evals[i];

	int counter = 0;
	for(i=0; i<N; i++) if(evals[i] == min_E) counter += 1;
	if(counter > 1) if(PRINT) printf("WARNING: GROUND STATE DEGENERACY!\n\n\n");


	delete[] evals, delete[] v_diag, delete[] ham_real;
	return min_E;
}



double cost(int N, double *state, double *hamiltonian){
	int INCX = 1,INCY = 1;
	double *state_conj, *state_temp, result;

	state_conj = new double[2*N]();
	state_temp = new double[2*N]();
	memcpy(state_conj, state, 2*N*sizeof(double)), memcpy(state_temp, state, 2*N*sizeof(double));


	if(CHECK) check_norm(state_temp, N), check_weights(state_temp, hamiltonian, N);

	matrix_vector_mult(hamiltonian, state_temp, N);//H*state, the operator acting on the ket, storing in state
	result = zdotc_(&N, state_conj, &INCX, state_temp, &INCY);//state* * state, the dot between the complex conj and the result of H*state

	delete[] state_conj, delete[] state_temp;
	return result;
}



int hop(unsigned long long int b, unsigned long long int *v,int n, int m){
	unsigned long long int i,x,y;
	int z_count = 0;

	x = (1ULL << (n-1)) + (1ULL << (m-1));
	for (i=n;i<m-1;i++)  if((1ULL<<i) & (b)) z_count++;
	y = (x ^ b);
	memcpy(v, &y, sizeof(unsigned long long int));
	return z_count%2;
}



int find(int N,unsigned long long int *v,unsigned long long int *b){
	int first = 0, last = N-1, mid;

	while (first <= last){
		mid = (int) ((first + last) / 2.0);
		if (*v > b[mid]) first = mid + 1;
		else if (*v < b[mid]) last = mid - 1;
		else return mid+1;
	}
}



unsigned long long int choose(int num_occupants){
	int i;
	unsigned long long int c=1ULL;
	for (i=0;i<num_occupants;i++) c=c*(NUMBER_OF_SITES-i);
	for (i=0;i<num_occupants;i++) c=c/(i+1);
	return c;
}



int combinations ( int N, int num_occupants,  unsigned long long int *b,int *tab){
	unsigned long long int x=0ULL,y;
	int i,c=0,d=0;

	for (i=0;i<num_occupants;i++) x=x+(1ULL<<i);
	b[0]=x;

	while ((c<NUMBER_OF_SITES)&& (d<num_occupants)){
		if (x & (1ULL<<c)) tab[d]=c+1, d++;
		c++;
	}

	for (i=1;i<N;i++){
		y = (x | (x - 1)) + 1;
		x = y | ((((y & -y) / (x & -x)) >> 1) - 1);
		b[i]=x;
		c=0, d=0;
		while ((c<NUMBER_OF_SITES)&& (d<num_occupants)){
			if (x & (1ULL<<c)) tab[i*num_occupants+d]=c+1, d++;
			c++;
		}
	}
}



double get_random_double(double lower, double upper, gsl_rng * rng){
	double zero_to_one, lower_to_upper;
	zero_to_one = gsl_rng_uniform(rng);
	lower_to_upper = zero_to_one*(upper-lower)+lower;
	return lower_to_upper;
}



double calc_distance(double initial, double target, double current){
	double distance;
	if(initial==target) return 0.0;
	distance = ((current - target) / (initial - target));
	return distance;
}


void calc_tau(Simulation_Parameters& sim_params){
	double tau_scalar;
	if(sim_params.tau == TAU_INIT){
		tau_scalar = 0.6/(1-sim_params.new_distance);
		if(tau_scalar > 2*TAU_SCALAR_BIG) tau_scalar = 2*TAU_SCALAR_BIG;
		if(tau_scalar < TAU_SCALAR) tau_scalar = TAU_SCALAR;
	}
	else if(sim_params.new_distance < 0.2) tau_scalar = TAU_SCALAR_TINY;
	else if((abs(sim_params.old_distance - sim_params.new_distance) < .1) and (sim_params.new_distance > 0.2 )) tau_scalar = TAU_SCALAR_BIG;
	else tau_scalar = TAU_SCALAR;

	sim_params.tau_previous = sim_params.tau;
	sim_params.tau = sim_params.tau*tau_scalar;
}



void construct_lattice(int lattice[NX][NY]){
	int x,y;
	for (x=0;x<NX;x++)
	{
		if (x%2 ==1) for (y=0;y<NY;y++) lattice[x][y] = (NX*x)+y+1;
		else for (y=0;y<NY;y++) lattice[x][NY-1-y] = NX*x+y+1;
	}
}



void assign_bonds(int *bonds, int lattice[NX][NY]){
	int bond_num=1, site,site2, i, neighbor_count, *neighbors;

	neighbors = new int[4]();

	for(site=1; site<NUMBER_OF_SITES+1;site++){
		neighbor_count = get_neighbors(site, neighbors, lattice);
		for(i=0;i<neighbor_count;i++){
			site2 = neighbors[i];
			if(site2>site){
				bonds[(site-1)*NUMBER_OF_SITES+site2-1]=bond_num;
				bonds[(site2-1)*NUMBER_OF_SITES+site-1]=bond_num;
				bond_num++;
			}
		}
	}
	delete[] neighbors;
}



int get_neighbors(int site, int *neighbors, int lattice[NX][NY]){
	int x,y,i,count=0;
	if (DEVICE_DIMENSION==2){
		for (x=0;x<NX;x++) for (y=0;y<NY;y++) if(site == lattice[x][y]) goto end_loop;//Finding the coordinates x and y
		end_loop: for(i=0;i<4;i++) neighbors[i]=0;

		if (PERIODIC){
			neighbors[0] = lattice[(x+1)%NX][y];
			neighbors[1] = lattice[(x+NX-1)%NX][y];
			neighbors[2] = lattice[x][(y+1)%NY];
			neighbors[3] = lattice[x][(y+(NY-1))%NY];
			count = 4;
		}else{
			if (x+1<NX) neighbors[count] = lattice[x+1][y], count++;
			if (x>0) neighbors[count] = lattice[x-1][y], count++;
			if (y+1<NY) neighbors[count] = lattice[x][y+1], count++;
			if (y>0) neighbors[count] = lattice[x][y-1], count++;
		}
		return count;
	}
	else if (DEVICE_DIMENSION==1){
		if (PERIODIC){
			neighbors[0] = (site%(NX*NY))+1;
			neighbors[1] = ((site+(NX*NY)-2)%(NX*NY))+1;
			count = 2;
		}else{
			if (site<(NX*NY)) neighbors[count] = site+1, count++;
			if (site>1) neighbors[count] = site-1, count++;
		}
		return count;
	}
	printf("\n\n\nERROR! UNSPECIFIED DIMENSION.\n\n\n");
	exit(0);
}


bool update_distances(Simulation_Parameters& sim_params){
	sim_params.temp_distance = sim_params.old_distance;
	sim_params.old_distance = sim_params.new_distance;
	sim_params.new_distance  = calc_distance(sim_params.initial_E, sim_params.ground_E,  sim_params.best_E);

	if(sim_params.new_distance > sim_params.old_distance && sim_params.sweeps_multiplier <= 2){
		if(PRINT){
			print_mc_results(sim_params);
			printf("############################################################################\n");
			printf("POOR CONVERGENCE DURING THE LAST SIMULATION FOR TAU: %f\nRUNNING AGAIN WITH TWICE THE SWEEPS\n", sim_params.tau);
			printf("############################################################################\n");
		}
		sim_params.new_distance  = 	sim_params.old_distance;
		sim_params.old_distance = sim_params.temp_distance;
		sim_params.sweeps_multiplier = sim_params.sweeps_multiplier * 2;
		return false;
	}

	sim_params.sweeps_multiplier = 1;
	return true;
}


void get_best_seed(Simulation_Parameters& sim_params){
	int INCX = 1,INCY = 1;

	for(sim_params.seed=1; sim_params.seed<NUM_SEEDS+1; sim_params.seed++){
		if(sim_params.E_array_fixed_tau[sim_params.seed-1] < sim_params.best_E){
			sim_params.best_E = sim_params.E_array_fixed_tau[sim_params.seed-1];
			memcpy(sim_params.best_evolved_state, &sim_params.evolved_state_fixed_tau[(sim_params.seed-1)*2*sim_params.N], 2*sim_params.N*sizeof(double));
			sim_params.evolved_target_dot_squared = pow(zdotc_(&sim_params.N, sim_params.target_state, &INCX, sim_params.best_evolved_state, &INCY),2);
		}
	}
}
