close all; clear all; clc 
% MOSEK test script for linear optimization
c = [3; 1; 5; 1];
a = [3 1 2 0;
     2 1 3 1;
     0 2 0 3];
blc = [30; 15; -inf];
buc = [30; inf; 25];
blx = zeros(4,1);
bux = [inf; 10; inf; inf];
[res] = msklpopt(c,a,blc,buc,blx,bux,[],'maximize');
sol = res.sol;
% Interior point solution 
x_sol = sol.itr.xx;
sol.itr.sux;
sol.itr.slx;
% Basic solution 
sol.bas.xx


%%
clear prob;

% Specify the c vector.
prob.c  = [3 1 5 1]';

% Specify a in sparse format.
% subi   = [1 1 1 2 2 2 2 3 3];
% subj   = [1 2 3 1 2 3 4 2 4];
% valij  = [3 1 2 2 1 3 1 2 3];
% 
% prob.a = sparse(subi,subj,valij);

prob.a = [3 1 2 0;
          2 1 3 1;
          0 2 0 3];
% Specify lower bounds of the constraints.
prob.blc = [30 15  -inf]';

% Specify  upper bounds of the constraints.
prob.buc = [30 inf 25 ]';

% Specify lower bounds of the variables.
prob.blx = zeros(4,1);

% Specify upper bounds of the variables.
prob.bux = [inf 10 inf inf]';

% Perform the optimization.
[r,res] = mosekopt('maximize',prob); 

% Show the optimal x solution.
res.sol.bas.xx
