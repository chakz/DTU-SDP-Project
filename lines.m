function [LDCincidenceMat, linesMatFrom, linesMatTo, nLines, linesFrom, linesTo, R, X, B, Z, lineMaxFlow, OriginBusLoc] = lines(mpc);
% Lines
nBuses = size(mpc.bus,1); % number of buses 
nLines = size(mpc.branch,1); % number of lines
%
OriginBusLoc = mpc.bus(:,1); % number of buses 
OriginLinesFrom = mpc.branch(:,1);
OriginLinesTo = mpc.branch(:,2);
[tf, linesFrom]=ismember(OriginLinesFrom,OriginBusLoc,'rows');
[tf, linesTo]=ismember(OriginLinesTo,OriginBusLoc,'rows');
%
%linesFrom = mpc.branch(:,1);
%linesTo = mpc.branch(:,2);
R = mpc.branch(:,3); % resistance
delta_R = 1e-4; 
R(R < delta_R) = delta_R; % adding small resistance to every tansformer with zero resistance
X = mpc.branch(:,4); % reactance
B = mpc.branch(:,5); % grosund admittance 
Z = R + j*X;
lineMaxFlow = mpc.branch(:,6); % rateA, (MVA long term rating)
% Line incidence matrix for SDP
linesMatFrom = sparse(1:nLines, linesFrom, 1, nLines, nBuses); % line incidence matrix [lines x nodes]
linesMatTo = sparse(1:nLines, linesTo, 1, nLines, nBuses); % line incidence matrix [lines x nodes]
% Line incidence matrix for LDC
LDCincidenceMat = zeros(nBuses,nLines); % bus-lines (node-arc) matrix with direction, element is 1 when the flow is from the bus and -1 when to the bus 
for i = 1:nLines
    lF = linesFrom(i);
    lT = linesTo(i);
    LDCincidenceMat(lF,i) = 1; 
    LDCincidenceMat(lT,i) = 1; 
end
end