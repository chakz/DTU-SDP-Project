function [YBus, Z] = busAdmittanceMatrix(R, j, X, B, nBuses, nLines, linesFrom, linesTo);
Z = R + j*X; % impedance [Ohms]
Y = 1./Z; % admittance [Siemens]
%B = j*B; 
YBus = zeros(nBuses,nBuses); % initiating the admittance matrix Ybus
%Formation of the off diagonal elements
for m = 1:nLines
    YBus(linesFrom(m),linesTo(m)) = YBus(linesFrom(m),linesTo(m)) - Y(m); % setting values for the right side
    YBus(linesTo(m),linesFrom(m)) = YBus(linesFrom(m),linesTo(m)); % mirroring the right side
end
%Formation of the diagonal elements
 for n = 1:nBuses
     for m = 1:nBuses
         if linesFrom(m) == n | linesTo(m) == n
             YBus(n,n) = YBus(n,n) + Y(m)+B(m);
         end
     end
 end
YBus;
end