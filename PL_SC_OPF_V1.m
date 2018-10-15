addpath(genpath('matpower6.0'))
mpc = IEEE_9BUS_Radial_modified;
%% AC OPF
mpopt = mpoption('model','AC', 'pf.tol', 1e-4,'opf.ac.solver','DEFAULT');
ACOPF = runopf(mpc,mpopt);
ACOPF_V = ACOPF.bus(:,8); % excluding slack bus (saving results for comparison) 
%% Defining Network Topology
[genMatrix,nGen, genLoc, baseMVA, PMin, PMax, QMin, QMax, nBuses, busLoc, Vmin, Vmax, Pd, Qd] = generators(mpc);
[LDCincidenceMat, linesMatFrom, linesMatTo, nLines, linesFrom, linesTo, R, X, B, Z, lineMaxFlow, OriginBusLoc] = lines(mpc);
% adds bus names if given in MatPower file 
if ~exist('mpc.bus_name','var')
    BusName = mpc.bus_name; 
    %BusInfo = [num2cell(busLoc'), num2cell(OriginBusLoc), BusName];
end
%% Bus Admittance Matrix
% [YBus, Z, susceptanceMatrix] = busAdmittanceMatrix(R, X, B, nBuses, nLines, linesFrom, linesTo); 
YBus_nplus1 = makeYbus(mpc); % with slack bus
YBus = YBus_nplus1(2:end,2:end); % without slack bus
ZBus_nplus1 = LDCincidenceMat*diag(Z);
ZBus = ZBus_nplus1(2:end,:);
%m = size(YBus,1); n = size(YBus,2); % YBus dimensions
%% Define parameters 
Vnom = ones(size(YBus,1),1); % nominal voltage vector [p.u.]
Vmin = ones(size(YBus,1),1)*0.95; % min voltage [p.u.]
Vmax = ones(size(YBus,1),1)*1.05; % max voltage [p.u.]
Sinj = ones(size(YBus,1),1)*0.1; % inverter capacity [kW]
Pav = ones(size(YBus,1),1)*0; % solar PV output [kW]
PF = 0.8; % power factor 
Pd = Pd(2:end);
Qd = Qd(2:end);
% zz = [1 1 1; 4 4 4; 7 7 7];
% zn = [1 2 3];
% zm = [1 2 3];
%%
%% Optimisation problem
cvx_begin 
    variable V(nBuses-1) complex; % Voltage vector [p.u.]
    variable Pc(nBuses-1); % Curtailed PV vector [kW]
    variable Pinj(nBuses-1); % Active power produced by inverter [kW]
    variable Qinj(nBuses-1); % Reative power absorbed/produced by inverter [kW] 
    Obj1 = 0;
    %for m = 1 : nLines
    %   for n = 1 : nLines
    %       Obj1 = Obj1 + (real(YBus(m,n))*real(V(m))+square(real(V(n)))); %Eq.1
    %   end
    %end
    Obj2 = 0;
    for n = 1 : size(YBus,1)
        %Obj2 = Obj2 + (80*Pc(n)+40*abs(Qinj(n))); % Eq.2
        Obj2 = Obj2 + (80*Pc(n)); % Eq.2
    end 
    Obj3 = 0; % Eq.3
    %for n = 1 : size(YBus,1)
    %    for l = 1 : size(YBus,1)
    %    Obj3 = Obj3 + (abs(V(n))-(1/((nBuses-1)+1))*sum(abs(V(l))))^2;
    %    end
    %end    
    minimize(Obj1+Obj2+Obj3) % Eq.4
       subject to 
       % Solar PV constraints
       for n = 1 : size(YBus,2)
            Pav(n) - Pinj(n) == Pc(n); % Eq.9
            0 <= Pinj(n) <= Pav(n); % Eq.10  
            square(Qinj(n)+Pinj(n)) <= Sinj(n)^2; % Eq.11
            Qinj(n) <= PF*Pinj(n); % Eq.12
       end  
       sub_realV = 0;
       sub_imagV = 0;
       sub_minV = 0;
       sub_maxV = 0;
       for n = 1 : size(YBus,2)
       % Voltage
           for l = 1 : size(YBus,1)
               sub_realV = sub_realV + real(ZBus(n,l))*(Pinj(l) - Pd(l)) + imag(ZBus(n,l))*(Qinj(l) - Qd(l)); % Eq.4
               sub_imagV = sub_imagV + imag(ZBus(n,l))*(Pinj(l) - Pd(l)) - real(ZBus(n,l))*(Qinj(l) - Qd(l)); % Eq.5
           end
           %real(V(n)) == Vnom(n) + sum(real(ZBus(n,:))*(Pinj(n) - Pd(n))+ imag(ZBus(n,:))*(Qinj(n) - Qd(n))); % Eq.5
           %imag(V(n)) == Vnom(n) + sum(imag(ZBus(n,:))*(Pinj(n) - Pd(n))- real(ZBus(n,:))*(Qinj(n) - Qd(n))); % Eq.6
           real(V(n)) == Vnom(n) + sub_realV; % Eq.4
           imag(V(n)) == Vnom(n) + sub_imagV; % Eq.5
       % Min voltage magnitude limit 
           for l = 1 : size(YBus,1)
               sub_minV = sub_minV + real(ZBus(n,l))*(Pinj(l) - Pd(l)) + imag(ZBus(n,l))*(Qinj(l) - Qd(l));
               sub_maxV = sub_maxV + real(ZBus(n,l))*(Pinj(l) - Pd(l)) + imag(ZBus(n,l))*(Qinj(l) - Qd(l)); 
           end
           Vmin(n) <= Vnom(n) + sub_minV; % Eq.7 
           Vnom(n) + sub_maxV <= Vmax(n);
       % Max voltage magnitude limit 
           %Vnom(n) + sum(real(ZBus(n,:))*(Pinj(n) - Pd(n))) + sum(imag(ZBus(n,:))*(Qinj(n) - Qd(n))) <= Vmax(n); % Eq.8  
       end  
cvx_end
% Results
Guggilam_V = [Vnom(1); real(V)]; % adding slack bus
VoltageTable = table(Guggilam_V, ACOPF_V, 'rownames', BusName) % compare Guggilam and ACOPF
VoltageTable.Properties.VariableNames = {'Guggilam' 'ACOPF'}
writetable(VoltageTable, 'V_Guggilam_VS_ACOPF', 'WriteRowNames',true)
% plot voltages 
plot(1:nBuses,[real(Vnom(1)); real(V)])
hold on
plot(1:nBuses,ACOPF_V)
xlabel('bus') 
ylabel('Voltage [p.u.]') 
xlim([1 nBuses]); ylim([0.9 1.1]) 
xticks(1:1:nBuses); 
legend('Guggilam', 'ACOPF')
set(gcf,'color','w'); 
