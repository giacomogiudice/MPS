# Matrix Product States 
by Giacomo Giudice

## Introduction
A set of functions for __MATLAB__ to work with _Matrix Product States_(MPS) in many-body quantum systems.
It is mostly focused on time evolution of systems using Suzuki-Trotter decomposition.

## Get Started
add a `addpath('<path>/mps') and you are good to go.

## GPU Support
If your GPU supports __CUDA__, you may try to speed up calculations by copying MPOs and MPSs to GPU, for example with `gpuArray(complex(<tensor>))`.
Functions `sweep`, `canonize`, `apply`, `contract` can work with objects distributed  on the GPU.
The other functions may not work, and one has to `gather` them back on the local workspace.  
Parallelism is beneficial mainly for large SVD decompositions, which is the main bottleneck.
Speedup will only be  obtained if the bond dimension, or the local Hilbert space are sufficiently large.  

