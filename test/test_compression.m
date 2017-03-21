% This script test the compression routines using both the SVD decomposition
% and the iterative compression from a random guess.
close all
clear all

%% Parameters
N = 18;
d = 3;

D_rand = 50;
D_force = 40;

truncation = 1e-6;
tolerance = 1e-8;

%% Compare Canonization Methods
Mrand = randomMPS(N,D_rand,d,1);
fprintf('Testing exact truncation...');
Msvd = sweep(Mrand,{},-1,D_rand+20,truncation);
[Miter,iterations] = sweep_iter(Mrand,{},randomMPS(N,D_rand+20,d,1),500,tolerance);
assert(1 - norm(braket(Mrand,Miter)) < tolerance);
assert(1 - norm(braket(Mrand,Msvd)) < 2*truncation);
fprintf('\tdone\n');
fprintf('\t distance original-SVD\t %.3g\n',1 - norm(braket(Mrand,Msvd)));
fprintf('\t distance original-iter\t %.3g\n',1 - norm(braket(Mrand,Miter)));
fprintf('\t distance SVD-iter %.3g\n',1 - norm(braket(Msvd,Miter)));



fprintf('Testing forced truncation...');
Msvd = sweep(Mrand,{},-1,D_force,truncation);
[Miter,iterations] = sweep_iter(Mrand,{},randomMPS(N,D_force,d,1),500,tolerance);
fprintf('\tdone\n');
fprintf('\t distance original-SVD %.3g\n',1 - norm(braket(Mrand,Msvd)));
fprintf('\t distance original-iter %.3g\n',1 - norm(braket(Mrand,Miter)));
fprintf('\t distance SVD-iter %.3g\n',1 - norm(braket(Msvd,Miter)));


fprintf('All tests passed!\n');