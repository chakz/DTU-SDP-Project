function [nLines, linesFrom, linesTo, R, X, B, j, SlineMax] = lines(mpc);
% Lines
nLines = size(mpc.branch,1); % number of lines
linesFrom = mpc.branch(:,1);
linesTo = mpc.branch(:,2);
R = mpc.branch(:,3); % resistance
X = mpc.branch(:,4); % reactance
B = mpc.branch(:,5); % ground admittance 
j = sqrt(-1); % imaginary part j
SlineMax = mpc.branch(:,6); % rateA, (MVA long term rating)
end