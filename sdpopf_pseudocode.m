addpath(genpath('matpower6.0'))
addpath(genpath('YALMIP-master')) 
%% Load the appropriate matpower case file
mpc = case9_SDP;
% runopf(mpc);
%% Defining Network Topology
[nLines, linesFrom, linesTo, R, X, B, j, SlineMax] = lines(mpc);
[nGen, genLoc, baseMVA, PMin, PMax, QMin, QMax, nBuses, busLoc, Vmin, Vmax, Pd, Qd] = generators(mpc);
% Network visualisation
network = digraph(linesFrom,linesTo); 
%network.Edges; % shows the number of Edges and Nodes 
%A = adjacency(network); % shows all nodes that llnes are connected to 
figure(1)
netgraph = plot(network,'Linewidth',2); 
highlight(netgraph,find(mpc.bus(:,3)>0),'NodeColor','red', 'MarkerSize',7); % loads 
highlight(netgraph,find(mpc.gen(:,1)>0),'NodeColor','black', 'MarkerSize',7); % gen
set(gcf,'color','w');
%% Generators incidence matrix
genMatrix = sparse(1:nGen, genLoc, 1 , nGen, nBuses); % gen incidence matrix  [plants x nodes]
figure(2)
spy(genMatrix,30) % plots the sparsity of matrix 
xlabel('bus')
ylabel('generator')
xlim([1 nBuses]); ylim([1 nGen])
xticks(1:1:nBuses); yticks(1:1:nGen)
set(gcf,'color','w');
% Lines incidence matrix
linesMatFrom = sparse(1:nLines, linesFrom, 1, nLines, nBuses); % line incidence matrix [lines x nodes]
linesMatTo = sparse(1:nLines, linesTo, 1, nLines, nBuses); % line incidence matrix [lines x nodes]
full(linesMatFrom)
%% Bus Admittance Matrix
[YBus,Z] = busAdmittanceMatrix(R, j, X, B, nBuses, nLines, linesFrom, linesTo);
%YBus = makeYbus(mpc);
%% Define your optimization variables
% Create the relevant vectors and assign them the appropriate size
% a vecor with number of lines. tranform into a matrix 
% e.g P = sdpvar(5, 1)
V = sdpvar(nBuses, 1,'full','complex'); % complex-valued fully parameterized vector
Vconj = sdpvar(nBuses, 1, 'full','complex'); 
W = sdpvar(2*nBuses, 2*nBuses,'hermitian','complex'); % hermitian matrix

%% Auxiliary variaibles 
ekMat = speye(nBuses); % sparse matrix with 1's on the main diagonal (e_k*e_k.')
Y_k =  ekMat*YBus; 
elmMat = adjacency(network); % sparse matrix with 1's corresponding to lines in/out
ellMat = speye(nBuses); % sparse matrix with 1's on the main diagonal 
Y_lm = ellMat*YBus - elmMat*YBus;

%%active power injection 
Yk = (1/2)*[real(Y_k + Y_k.') imag(Y_k.' - Y_k); 
            imag(Y_k - Y_k.') real(Y_k + Y_k.')]; 
% reactive power injection
Yk_ = -(1/2)*[imag(Y_k + Y_k.') real(Y_k - Y_k.'); 
              real(Y_k.' - Y_k) imag(Y_k + Y_k.')]; 
% matrix for the square of voltage magnitude          
Mk = blkdiag(ekMat, ekMat); 
% line active power flow
Ylm = (1/2)*[real(Y_lm + Y_lm.') imag(Y_lm.' - Y_lm); 
            imag(Y_lm - Y_lm.') real(Y_lm + Y_lm.')]; 
% line reactive power injection 
Ylm_ = -(1/2)*[real(Y_lm + Y_lm.') imag(Y_lm.' - Y_lm); 
            imag(Y_lm - Y_lm.') real(Y_lm + Y_lm.')]; 
%% Extension (D.Molzahn 2013) - makesdpmat (YALMIP)
%function [Yk,Yk_,Mk,Ylineft,Ylinetf,Y_lineft,Y_linetf,YL,YL_] = makesdpmat(mpc)
% allowing multiple generators on the same bus
% allowing parallel lines 
% 
% lines settings
% ratio = mpc.branch(:,9); % ratio, transformer off nominal turns ratio
% ratio(ratio == 0) = 1; % avoiding division by zero
% theta = mpc.branch(:,10)*pi/180; % angle, transformer phase shift angle (radians)
% 
% gbcosft = real(Z).*cos(theta) + imag(Z).*cos(theta+pi/2);
% gbsinft = real(Z).*sin(theta) + imag(Z).*sin(theta+pi/2);
% gbcostf = real(Z).*cos(-theta) + imag(Z).*cos(-theta+pi/2);
% gbsintf = real(Z).*sin(-theta) + imag(Z).*sin(-theta+pi/2);

% Active power flow from the original bus to another bus
% YlineFT_P = 0.5*(sparse(  [linesFrom linesFrom linesFrom linesFrom+nBuses linesFrom+nBuses linesFrom+nBuses], ...         % from
%                           [linesFrom linesTo linesTo+nBuses linesFrom+nBuses linesTo linesTo+nBuses], ...                 % to 
%                           [real(Z)./(ratio.^2) -gbcosft./ratio gbsinft./ratio real(Z)./(ratio.^2) -gbsinft./ratio -gbcosft./ratio] ...% 
%                            ,2*nBuses,2*nBuses) + ... 
%                sparse(    [linesFrom linesFrom linesFrom linesFrom+nBuses linesFrom+nBuses linesFrom+nBuses], ...
%                           [linesFrom linesTo linesTo+nBuses linesFrom+nBuses linesTo linesTo+nBuses], ...
%                           [real(Z)./(ratio.^2) -gbcosft./ratio gbsinft./ratio real(Z)./(ratio.^2) -gbsinft./ratio -gbcosft./ratio] ...
%                            ,2*nBuses,2*nBuses).'); % number of busses
% % Reactive power flow from another bus back to the original bus                       
% YlineFT_Q = 0.5*(sparse(   [linesFrom linesFrom linesFrom linesFrom+nBuses linesFrom+nBuses linesFrom+nBuses], ...
%                            [linesFrom linesTo linesTo+nBuses linesFrom+nBuses linesTo linesTo+nBuses], ...
%                            [-(imag(Z)+B./2)./(ratio.^2) gbsinft./ratio gbcosft./ratio -(imag(Z)+B./2)./(ratio.^2) -gbcosft./ratio gbsinft./ratio] ...
%                             ,2*nBuses,2*nBuses) + ... 
%                 sparse(    [linesFrom linesFrom linesFrom linesFrom+nBuses linesFrom+nBuses linesFrom+nBuses], ...
%                            [linesFrom linesTo linesTo+nBuses linesFrom+nBuses linesTo linesTo+nBuses], ...
%                            [-(imag(Z)+B./2)./(ratio.^2) gbsinft./ratio gbcosft./ratio -(imag(Z)+B./2)./(ratio.^2) -gbcosft./ratio gbsinft./ratio] ...
%                             ,2*nBuses,2*nBuses).');
% % Active power flow from another bus back to the original bus
% YlineTF_P = 0.5*(sparse(  [linesFrom linesFrom linesFrom+nBuses linesFrom+nBuses linesTo linesTo+nBuses], ...
%                           [linesTo linesTo+nBuses linesTo linesTo+nBuses linesTo linesTo+nBuses], ...
%                           [-gbcostf./ratio -gbsintf./ratio gbsintf./ratio -gbcostf./ratio real(Z) real(Z)] ...
%                            ,2*nBuses,2*nBuses) + ... 
%                sparse(    [linesFrom linesFrom linesFrom+nBuses linesFrom+nBuses linesTo linesTo+nBuses], ...
%                           [linesTo linesTo+nBuses linesTo linesTo+nBuses linesTo linesTo+nBuses], ...
%                           [-gbcostf./ratio -gbsintf./ratio gbsintf./ratio -gbcostf./ratio real(Z) real(Z)] ...
%                            ,2*nBuses,2*nBuses).');
% % Reactive power flow from another bus back to the original bus                        
% YlineTF_Q = 0.5*(sparse(   [linesFrom linesFrom linesFrom+nBuses linesFrom+nBuses linesTo linesTo+nBuses], ...
%                            [linesTo linesTo+nBuses linesTo linesTo+nBuses linesTo  linesTo+nBuses], ...
%                            [gbsintf./ratio -gbcostf./ratio gbcostf./ratio gbsintf./ratio -(imag(Z)+B./2) -(imag(Z)+B./2)] ...
%                             ,2*nBuses,2*nBuses) + ...
%                 sparse(    [linesFrom linesFrom linesFrom+nBuses linesFrom+nBuses linesTo linesTo+nBuses], ...
%                            [linesTo linesTo+nBuses linesTo linesTo+nBuses linesTo linesTo+nBuses], ...
%                            [gbsintf./ratio -gbcostf./ratio gbcostf./ratio gbsintf./ratio -(imag(Z)+B./2) -(imag(Z)+B./2)] ...
%                             ,2*nBuses,2*nBuses).');
%                         
% Matrices to calculate active and reactive power line losses
% YL = sparse([linesFrom linesFrom linesFrom+nBuses linesFrom+nBuses linesTo linesTo linesTo+nBuses linesTo+nBuses], ...
%             [linesFrom linesTo linesFrom+nBuses linesTo+nBuses linesFrom linesTo linesFrom+nBuses linesTo+nBuses], ...
%             [ones(nBuses,1) -ones(nBuses,1) ones(nBuses,1) -ones(nBuses,1) ones(nBuses,1) -ones(nBuses,1) ones(nBuses,1) -ones(nBuses,1)] ...
%              ,2*nBuses,2*nBuses)*R* (real(Z).^2+imag(Z).^2);
%          
% YL_ = sparse([linesFrom  linesFrom linesFrom+nBuses linesFrom+nBuses linesTo linesTo linesTo+nBuses linesTo+nBuses], ...
%              [linesFrom linesTo linesFrom+nBuses linesTo+nBuses linesTo linesTo linesFrom+nBuses linesTo+nBuses], ...
%              [ones(nBuses,1) -ones(nBuses,1) ones(nBuses,1) -ones(nBuses,1) -ones(nBuses,1) ones(nBuses,1) -ones(nBuses,1) ones(nBuses,1)] ...
%               ,2*nBuses,2*nBuses)* X.* (real(Z).^2+imag(Z).^2) + ...
%      -sparse([linesFrom linesFrom+nBuses linesTo linesTo+nBuses], ...
%              [linesFrom linesFrom+nBuses linesTo linesTo+nBuses], ...
%              [ones(nBuses,1) ones(nBuses,1) ones(nBuses,1) ones(nBuses,1)] ...
%               ,2*nBuses,2*nBuses) * B./2;
%% Objective function
ck0 = mpc.gencost(:,7); % gen cost coefficient
ck1 = mpc.gencost(:,6)*baseMVA;
ck2 = mpc.gencost(:,5)*baseMVA^2;
Objective = sum(ck2.*trace(Yk.*W).^2 + ck1.*trace(Yk.*W) + ck0); % Objective = c2*P^2+c1*P+c0


%% Constraints 
% Hint: To improve numerical performance of the solver, define inequality
% constraints with identical lower and upper bounds as equality constraints
% see https://yalmip.github.io/tutorial/basics/

%Constraints = [trace(Yk*W) (P) <= Pmax - Pd,...
%               trace(Yk*W) (P) >= Pmin - Pd,...
%               trace(Yk_*W) (Q) <= Qmax - Qd,...
%               trace(Yk_*W) (Q) >= Qmin - Qd,...
%               trace(Mk*W) (V) >= Vmax.^2,... % need Schur's complement 
%               trace(Mk*W) (V) >= Vmin.^2,... % need Schur's complement 
%               trace(Ylm*W)^2+trace(Ylm_*W)^2 <= SlineMax^2,...
        %inequality constraints
        %upper and lower bounds on optimization variables
                %]
            
% 1a-1f from Lavaei paper! 
            
%% Run the optimization
%optimize(Constraints, Objective, sdpsettings('solver','gurobi'))

%% Print the results
% e.g. value(P)


% Calculate active and reactive power injections


% Compare to non-convex AC-OPF results


%% Post-process
% eig(A) might be helpful to calculate eigenvalues and -vectors of matrix A


% Calculate eigenvalue ratio and evaluate exactness of the relaxation


% Decompose W matrix to obtain optimal voltage vector

%% Vary network properties and evaluate exactness of relaxation



%% 3 bus test case: Load the appropriate matpower case file

% e.g. mpc = case9_SDP

%% Solve optimization problem


%% Change line limits


%% Bonus: Include penalty term in objective function




