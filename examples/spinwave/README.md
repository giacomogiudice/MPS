## Example: Spin Waves in the Heisenberg Model

This example demonstrates time-evolution through a variety of different methods.
In particular it shows the propagation of spin waves in a one-dimensional lattice with antiferromagnetic nearest neighbor interactions, i.e. the Heisenberg model described with the Pauli matrices  
_H = ∑<sub>k</sub> σ<sup>x</sup><sub>k</sub>σ<sup>x</sup><sub>k+1</sub> + σ<sup>y</sup><sub>k</sub>σ<sup>y</sup><sub>k+1</sub> + σ<sup>z</sup><sub>k</sub>σ<sup>z</sup><sub>k+1</sub>_  
The initial state is chosen to be _|100...000⟩_.


The parameters of the simulation are defined in `setup.m`. To run some simulations start by adding the path from your current working directory. For example from the trunk of this project run from your __MATLAB__ prompt  
`addpath(genpath('examples'))` or `addpath('examples/spinwave')`   
for this specific directory. Run a simulation with   
`spinwave_<name>`  
The results should be saved by default in the `workspaces` folder. To compare the results run  
`plot_comparison`  
You should see an animation of the magnetization profiles over time.
