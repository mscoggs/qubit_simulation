# qubit_sim
Simulating the evolution on a superconducting chip, with the goal of finding patterns in the protocols which can inform ansatz for largers system sizes. In data_generation there are 4 different types of simulation; adiabatic, monte-carlo brute-force (mcbf), monte-carlo bang-bang (mcbb), and monte-carlo discrete-bang (mcdb). The first uses a linear parametrization for the protocols. The seconds uses a discrete approximation, allowing the protocols to take on any value on a given interval. The third assumes bang-bang protocols, and the parameters becomes the times that transtions occur (going from 1 to 0 or 0 to 1). The fourth is like the mcbf, except the values are restricted to 1 or 0. 

Set up the simulation parameters in parameters.h. Here you specify which simulation(s) you want to run along with setting some of the simulation parameters. Run 'make compile', then ./main occupancy_number initial_j initial_k target_j target_k to generate the data.

## REQUIREMENTS FOR DATA GENERATION
- c++
- g++ compiler
- blas
- lapack
- gsl


## REQUIREMENTS FOR DATA ANALYSIS
 - anaconda
     -jupyter notebook
     -python
     -glob 
     -scipy 
     -matplotlib 
     -mplot3d 
     -pandas 
 
## PAPER, UNDER REVIEW
[Topological and geometric patterns in optimal bang-bang protocols for variational quantum algorithms: application to the XXZ model on the square lattice](https://arxiv.org/abs/2012.05476)
