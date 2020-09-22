#include "check.h"
#include "hamiltonian.h"
#include "operations.h"
#include "parameters.h"



void construct_device_hamiltonian(Simulation_Parameters sim_params, double *hamiltonian ,double *j_array, double *k_array, double *b_array,  int index){
	int i,ii,j,state,site,sign,bond,neighbor_count,*neighbors;
	unsigned long long int *v,comparison=0;

	v         = new unsigned long long int[1]();
	neighbors = new int[4]();

	for (i=0; i<sim_params.N*sim_params.N*2; i++) hamiltonian[i] = 0.0;

	for (i=0;i<sim_params.N;i++){

		for (j=0;j<sim_params.num_occupants;j++){//The J term calculation
			site=sim_params.table[i*sim_params.num_occupants+j];
			neighbor_count = get_neighbors(site, neighbors, sim_params.lattice);
			for (ii=0; ii<neighbor_count; ii++){
				if (((1ULL<<(neighbors[ii]-1))&sim_params.b[i])==0){//making sure neighbor is not occupied, otherwise nothing happens

					hop(sim_params.b[i], v,site, neighbors[ii]);
					state=find(sim_params.N,v, sim_params.b);
					bond = sim_params.bonds[(site-1)*NUMBER_OF_SITES+neighbors[ii]-1]-1;
					hamiltonian[(sim_params.N*i+state-1)*2] -= j_array[index*NUMBER_OF_SITES*2+bond];
				}
			}
		}
		for (j=1;j<(NUMBER_OF_SITES);j++){//The K term calculation
			site=j;
			neighbor_count = get_neighbors(site, neighbors, sim_params.lattice);
			for (ii=0; ii<neighbor_count;ii++){
				if (neighbors[ii] > site){
					sign = -1;
					comparison = (1ULL<<(neighbors[ii]-1))+(1ULL<<(site-1));
					if((comparison&sim_params.b[i])==comparison || (comparison&sim_params.b[i])==0) sign = 1;
					bond = sim_params.bonds[(site-1)*NUMBER_OF_SITES+neighbors[ii]-1]-1;
					hamiltonian[((sim_params.N*i)+i)*2] += k_array[index*NUMBER_OF_SITES*2+bond]*sign;
				}
			}
		}
		/*for (j=0; j<NUMBER_OF_SITES;j++){//The sim_params.b term calculation
			sign = -1;
			if(((1ULL<<j)&sim_params.b[i])>0) sign=1;
			hamiltonian[((sim_params.N*i)+i)*2] += b_array[index*NUMBER_OF_SITES+j]*sign;
		}*/
	}
	delete[] neighbors, delete[] v;
	if(CHECK) check_hermicity(hamiltonian, sim_params.N);
}



void construct_device_hamiltonian_uniform(Simulation_Parameters sim_params, double *hamiltonian ,double *jkb){
	int i,ii,j,state,site,sign,neighbor_count,*neighbors;
	unsigned long long int *v,comparison=0;

	v         = new unsigned long long int[1]();
	neighbors = new int[4]();

	for (i=0; i<sim_params.N*sim_params.N*2; i++) hamiltonian[i] = 0.0;

	for (i=0;i<sim_params.N;i++){

		for (j=0;j<sim_params.num_occupants;j++){//The J term calculation

			site=sim_params.table[i*sim_params.num_occupants+j];
			neighbor_count = get_neighbors(site, neighbors, sim_params.lattice);

			for (ii=0; ii<neighbor_count; ii++){
				if (((1ULL<<(neighbors[ii]-1))&sim_params.b[i])==0){//making sure neighbor is not occupied, otherwise nothing happens
					hop(sim_params.b[i], v,site, neighbors[ii]);
					state=find(sim_params.N,v, sim_params.b);
					hamiltonian[(sim_params.N*i+state-1)*2] -= jkb[0];
				}
			}
		}
		for (j=1;j<(NUMBER_OF_SITES);j++){//The K term calculation
			site=j;
			neighbor_count = get_neighbors(site, neighbors, sim_params.lattice);
			for (ii=0; ii<neighbor_count;ii++){
				if (neighbors[ii] > site){
					sign = -1;
					comparison = (1ULL<<(neighbors[ii]-1))+(1ULL<<(site-1));
					if((comparison&sim_params.b[i])==comparison || (comparison&sim_params.b[i])==0) sign = 1;
					hamiltonian[((sim_params.N*i)+i)*2] += jkb[1]*sign;
				}
			}
		}
		/*for (j=0; j<NUMBER_OF_SITES;j++){//The B term calculation
			sign = -1;
			if(((1ULL<<j)&sim_params.b[i])>0) sign=1;
			hamiltonian[((sim_params.N*i)+i)*2] += jkb[2]*sign;
		}*/
	}
	delete[] neighbors, delete[] v;
	if(CHECK) check_hermicity(hamiltonian, sim_params.N);
}



void construct_model_hamiltonian(Simulation_Parameters sim_params, double *hamiltonian){
	int i,ii,j,state,site,neighbor,sign,neighbor_count, *neighbors ;
	unsigned long long int *v,comparison;
	int T=1, V=1;
	v         = new unsigned long long int[1]();
	neighbors = new int[4]();

	for (i=0; i<sim_params.N*sim_params.N*2; i++) hamiltonian[i] = 0;

	for (i=0;i<sim_params.N;i++){
		for (j=0;j<sim_params.num_occupants;j++){//The T term calculation

			site=sim_params.table[i*sim_params.num_occupants+j];
			neighbor_count = get_neighbors(site, neighbors, sim_params.lattice);

			for (ii=0; ii<neighbor_count; ii++){
				if (((1ULL<<(neighbors[ii]-1))&sim_params.b[i])==0){//checking if the neighbor is occupied

					sign = hop(sim_params.b[i], v,site, neighbors[ii]);
					if (sign==0) sign=1;
					else sign=-1;
					state=find(sim_params.N,v, sim_params.b);
					hamiltonian[(sim_params.N*i+state-1)*2] -= (T*sign);
				}
			}
		}
		for (j=1;j<(NUMBER_OF_SITES);j++){//The V term calculation
			site=j;
			neighbor_count = get_neighbors(site, neighbors, sim_params.lattice);
			for (ii=0; ii<neighbor_count;ii++){
				if (neighbors[ii] > site){
					sign = -1;
					comparison = (1ULL<<(neighbors[ii]-1))+(1ULL<<(site-1));
					if((comparison&sim_params.b[i])==comparison || (comparison&sim_params.b[i])==0) sign = 1;
					hamiltonian[((sim_params.N*i)+i)*2] += sign*V;
				}
			}
		}
	}
	delete[] neighbors, delete[] v;
	if(CHECK) check_hermicity(hamiltonian, sim_params.N);
}
