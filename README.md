# Matrix Product States 
by Giacomo Giudice

## Notes
This has been tested using __MATLAB__ R2016a down to R2014a. 
Given Mathwork's history of releasing non-backwards-compatible updates, things will probably break in future versions.  

## Introduction
A (very) small library for to simulate many-body quantum systems with _Matrix Product States_ (MPS).
It is designed to be scalable and performant at the same time staying flexible enough to be hackable.
It is mostly focused on time evolution of _finite_ MPS, both through Time-Evolving Block Decimation (TEBD) using the Trotter-Suzuki decomposition, as well as the Time-Dependent Variational Principle (TDVP).
It also features single-site Iterative Variational Optimization to find ground states of Hamiltonian systems.

## Features
* __Scalable__ All routines are written with performance in mind, all underlying contractions are optimal for large bond dimension _D_. The complexity should never be more than _O(D<sup>3</sup>)_.
* __Efficient__ The underlying and most time-consuming operations, such as SVD decomposition and tensor contraction reduce to built-in __MATLAB__ operations, which are very fast.
* __Hackable__ The different matrix-product objects are not hidden under many layers of abstraction, but are simple cell arrays. Objects can then be easily inspected, and the code is meant to be played around with.

## Getting Started
Add a `addpath('<path>/mps')` and you are good to go.
A good place to start is the `examples` folder. There are currently some MWEs of different methods for time-evolution in `examples/spinwave` and some examples of ground state estimation routines in `examples/ground_states`.
Additionally, the `test` folder demonstrates the use of the elementary functions.

## GPU Support
__CUDA__ support has been discontinued. Copying matrix product elements to GPU is possible with `gpuArray(complex(<tensor>))`, but correct memory management is critical.

## References
1. U. Schollwöck, _The density-matrix renormalization group in the age of matrix product states_, Annals of Physics (2011).
2. F. Verstraete, V. Murg, and J.I. Cirac, _Matrix product states, projected entangled pair states, and variational renormalization group methods for quantum spin systems_, Advances in Physics (2008).
3. J. Haegeman, C. Lubich, I. Oseledets, B. Vandereycken, F. Verstraete, _Unifying time evolution and optimization with matrix product states_, Physical Review B (2016).
