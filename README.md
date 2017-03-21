# Matrix Product States 
by Giacomo Giudice

## Introduction
A set of functions for __MATLAB__ to work with _Matrix Product States_ (MPS) in many-body quantum systems.
It is mostly focused on time evolution (TEBD) using Suzuki-Trotter decomposition. 

## Getting Started
Add a `addpath('<path>/mps')` and you are good to go.
The `test` folder is probably a good place to start looking at examples.


## GPU Support

__CUDA__ support has been discontinued. Copying matrix product elements to GPU is possible with `gpuArray(complex(<tensor>))`, but correct memory management is critical.

## References
1. Schollwöck, U., 2011. The density-matrix renormalization group in the age of matrix product states. _Annals of Physics_, 326(1), pp.96-192.
2. Verstraete, F., Murg, V. and Cirac, J.I., 2008. Matrix product states, projected entangled pair states, and variational renormalization group methods for quantum spin systems. _Advances in Physics_, 57(2), pp.143-224.
