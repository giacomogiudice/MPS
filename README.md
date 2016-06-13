# Matrix Product States 
by Giacomo Giudice

## Introduction
A set of functions for __MATLAB__ to work with _Matrix Product States_ (MPS) in many-body quantum systems.
It is mostly focused on time evolution of systems using Suzuki-Trotter decomposition.

## Getting Started
Add a `addpath('<path>/mps')` and you are good to go.
The `test` folder is probably a good place to start looking at examples.

## Routines
In the `routines` folder one can find different methods to simulate driven-dissipative Bose--Hubbard chains.
Those include:
* Lindbladian Evolution
* Lindbladian Inversion
* Integration of the Stochastic Schrödinger Equation
* Evolution of the density matrix by encoding it as an MPS.

## GPU Support
If your GPU supports __CUDA__, you may try to speed up calculations by copying MPOs and MPSs to GPU, for example with `gpuArray(complex(<tensor>))`.
Functions `sweep`, `canonize`, `apply`, `contract` can work with objects distributed  on the GPU.
The other functions may not work, and one has to `gather` them back on the local workspace.  
Parallelism is beneficial mainly for large SVD decompositions, which is the main bottleneck.
Speedup will only be obtained if the product of the bond dimension and the local Hilbert space dimension are sufficiently large.  
See `routines/density_mps.m` and `util/dm_gpu.m' too see an example.  

## References
1. Schollwöck, U., 2011. The density-matrix renormalization group in the age of matrix product states. _Annals of Physics_, 326(1), pp.96-192.
2. Verstraete, F., Murg, V. and Cirac, J.I., 2008. Matrix product states, projected entangled pair states, and variational renormalization group methods for quantum spin systems. _Advances in Physics_, 57(2), pp.143-224.
