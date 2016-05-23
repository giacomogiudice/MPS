# Matrix Product States 
by Giacomo Giudice

## Introduction
A set of functions for __MATLAB__ to work with _Matrix Product States_ (MPS) in many-body quantum systems.
It is mostly focused on time evolution of systems using Suzuki-Trotter decomposition.

## Get Started
Add a `addpath('<path>/mps')` and you are good to go.

## GPU Support
If your GPU supports __CUDA__, you may try to speed up calculations by copying MPOs and MPSs to GPU, for example with `gpuArray(complex(<tensor>))`.
Functions `sweep`, `canonize`, `apply`, `contract` can work with objects distributed  on the GPU.
The other functions may not work, and one has to `gather` them back on the local workspace.  
Parallelism is beneficial mainly for large SVD decompositions, which is the main bottleneck.
Speedup will only be  obtained if the bond dimension, or the local Hilbert space are sufficiently large.  

## References
1. U. Schollwöck. The density-matrix renormalization group in the age of matrix product states. _Annals of Physics_ 326.1 (2011): 96-192.  
2. F. Verstraete, J. J. Garcia-Ripoll and J. I. Cirac. Matrix product density operators: simulation of finite-temperature and dissipative systems. _Physical Review Letters_ 93.20 (2004): 207204.
