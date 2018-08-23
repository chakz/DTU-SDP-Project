function [nGen, genLoc, baseMVA, PMin, PMax, QMin, QMax, nBuses, busLoc, Vmin, Vmax, Pd, Qd] = generators(mpc);
nGen = size(mpc.gen,1); % number of generators
genLoc = mpc.gen(:,1); % location of generators
baseMVA = mpc.baseMVA; % power rating
PMin = mpc.gen(:,10)/baseMVA; % 
PMax = mpc.gen(:,9)/baseMVA;
QMin = mpc.gen(:,4)/baseMVA;
QMax = mpc.gen(:,5)/baseMVA;
% Buses
nBuses = size(mpc.bus,1); % number of buses 
busLoc = mpc.bus(:,1); % number of buses 
Vmin = mpc.bus(:,13);
Vmax = mpc.bus(:,12);
% Demand 
Pd = mpc.bus(:,3)/baseMVA; % active power
Qd = mpc.bus(:,4)/baseMVA; % reactive power
end