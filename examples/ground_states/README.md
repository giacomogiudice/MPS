## Example: Variational Ground State Search

In this example we find the ground state of some different models. In each routine, the Hamiltonian is constructed as an MPO, and the ground state is found by applying the variational principle iteratively. Starting from a random MPS ansatz with a fixed bond dimension, the problem is reduced to a local optimization problem. 
This example demonstrates time-evolution through a variety of different methods.
From the trunk add the correct path  
`addpath(genpath('examples'))` 
and then run a script inside it, such as  
`heisenberg_gs`  

The different models available are  
* __the Heisenberg model__ with an external field described by _H = ∑<sub>k</sub> σ<sup>x</sup><sub>k</sub>σ<sup>x</sup><sub>k+1</sub> + σ<sup>y</sup><sub>k</sub>σ<sup>y</sup><sub>k+1</sub> + σ<sup>z</sup><sub>k</sub>σ<sup>z</sup><sub>k+1</sub> + gσ<sup>x</sup><sub>k</sub>_, where _σ<sup>a</sup>_ are the Pauli matrices. The simulation can be run with `heisenberg_gs`. At the end of the routine, the magnetization along _z_ and the _z-z_ correlation are plotted.
* __the Bose-Hubbard model__ described by _H = ∑<sub>k</sub> -µa<sup>†</sup><sub>k</sub>a<sub>k</sub> + a<sup>†</sup><sub>k</sub>a<sup>†</sup><sub>k</sub>a<sub>a<sub>k</sub> - Ja<sup>†</sup><sub>k</sub>a<sub>a<sub>k+1</sub> - -Ja<sub>k</sub>a<sub>a<sup>†</sup><sub>k+1</sub>_. The simulation can be run with `bose_hubbard_gs`. The occupation number and the _g<sup>2</sup>_ correlation are plotted.

