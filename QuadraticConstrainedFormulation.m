close all; clear all; clc 
% Quadratic Objective formulation 
% minimize x1^2 + 0.1x2^2 + x3^2 -x1x3 - x2
% s.t.c 1 <= x1 + x2 + x3
%       0 <= x

q = [2 0 -1;
     0 0.2 0;
     -1 0 2];
 
c = [0; -1; 0];

a = [1 1 1];

blc = 1.0;

buc = inf;

%blx = zeros(3,1);
blx = sparse(3,1);
%bux = [inf; inf; inf];
bux = [];
% Optimization function 
res = mskqpopt(q, c, a, blc, buc, blx, bux);

% The results are as follows
res.sol.itr.xx

% Quadratic Constraint Formulation
close all; clear all; clc 
clear prob;

% Specify the linear objective terms.
prob.c      = [0, -1, 0];

% Specify the quadratic terms of the constraints.
prob.qcsubk = [1     1    1   1  ]';
prob.qcsubi = [1     2    3   3  ]';
prob.qcsubj = [1     2    3   1  ]';
prob.qcval  = [-2.0 -2.0 -0.2 0.2]';

% Specify the quadratic terms of the objective.
prob.qosubi = [1     2    3    3  ]';
prob.qosubj = [1     2    3    1  ]';
prob.qoval  = [2.0   0.2  2.0 -1.0]';

% Specify the linear constraint matrix
prob.a      = [1 1 1];

% Specify the lower bounds
prob.blc    = [1];
prob.blx    = zeros(3,1);

[r,res]     = mosekopt('minimize',prob);

% Display the solution.
fprintf('\nx:');
fprintf(' %-.4e',res.sol.itr.xx');
fprintf('\n||x||: %-.4e',norm(res.sol.itr.xx));