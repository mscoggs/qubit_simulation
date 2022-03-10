# qubit_sim
Simulating the evolution of a superconducting chip with the goal of finding patterns in the optimal protocols (values of the controls over time which evolve an initial state into a target state in the shortest possible time) over a variety of initial and target combinations. These patterns can inform ansatz for largers system sizes, making them computationally feasible or allowing the use of a hybrid quantum device where evolution is carried out on a chip. We focus on the XXZ model.

In data_generation there are 4 different types of simulation which all try to find the smallest total time required, given their parameterizations, for evolving an initial state into a target state:
 - Adiabatic (ADIA): evolves that state according to a linear parameterization, incrementing time until the target state is achieved.
 - Monte-Carlo Brute-Force (MCBF): This breaks up the protocol into N even intervals, allowing the protocol to take on any value on a given interval. This is very inefficient but is useful for confirming the bang-bang nature of the protocols.
 - Monte-Carlo Bang-Bang (MCBB): We assume bang-bang protocols (where the control is at its maximum or minimum at any given time) and make the parameter the time that the transitions occur. This is much more efficient the the MCBF and much more accurate.
 - Monte-Carlo Discrete-Bang (MCDB): This is similar to the MCBF except the intervals are restricted to either 1 or 0. This is nearly as accurate as the MCBB but much more efficient. There is an option to run a MCBB secondary for each total time, which is the optimal combination of accuracy and computational efficiency -- all of our available data was generated using this.



# data_generation

## REQUIREMENTS
- c++
- g++ compiler
- blas
- lapack
- gsl

## USAGE
Set up the simulation parameters in parameters.h. Here you specify which simulation(s) you want to run along with setting some of the simulation parameters. After that:
```bash  
make compile
./main occupancy_number initial_j initial_k target_j target_k
```
generates the data for this initial-target combination with hamiltonian parameters j and k.



# data_analysis

## REQUIREMENTS
 - anaconda
     - jupyter notebook
     - python
     - glob 
     - scipy 
     - matplotlib 
     - mplot3d 
     - pandas 
 
 ## USAGE
 ```bash  
jupyter notebook analyze_data.ipynb
```
then follow the markdown in the notebook



# Paper
[Topological and geometric patterns in optimal bang-bang protocols for variational quantum algorithms: application to the XXZ model on the square lattice](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.3.043165)
