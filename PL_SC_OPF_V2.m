clear; clc;
addpath(genpath('matpower6.0'))
addpath(genpath('matpower Project'))
mpc = IEEE_9BUS_Radial_der_modified;
%mpc = case33bw;
%[new_mpc, idx] = addgen2mpc(mpc, 3, 1, 'solar')
%% AC OPF
mpopt = mpoption('model','AC', 'pf.tol', 1e-4,'opf.ac.solver','DEFAULT');
ACOPF = runopf(mpc,mpopt);
ACOPF_V = ACOPF.bus(:,8); % excluding slack bus (saving results for comparison) 
%% Defining Network Topology
[genMatrix,nGen, genLoc, baseMVA, PMin, PMax, QMin, QMax, nBuses, busLoc, Vmin, Vmax, Pd, Qd] = generators(mpc);
[LDCincidenceMat, linesMatFrom, linesMatTo, nLines, linesFrom, linesTo, R, X, B, Z, lineMaxFlow, OriginBusLoc] = lines(mpc);
% adds bus names if given in MatPower file 
if exist('mpc.bus_name','var') == 1
    BusName = mpc.bus_name; 
    %BusInfo = [num2cell(busLoc'), num2cell(OriginBusLoc), BusName];
end
%% Bus Admittance Matrix
% [YBus, Z, susceptanceMatrix] = busAdmittanceMatrix(R, X, B, nBuses, nLines, linesFrom, linesTo); 
YBus_nplus1 = makeYbus(mpc); % with slack bus
YBus = YBus_nplus1(2:end,2:end); % without slack bus
ZBus_nplus1 = LDCincidenceMat*diag(Z);
ZBus = ZBus_nplus1(2:end,:);
%% Define parameters 
Vnom = ones(size(YBus,1),1); % nominal voltage vector [p.u.]
Vmin = ones(size(YBus,1),1)*0.95; % min voltage [p.u.]
Vmax = ones(size(YBus,1),1)*1.05; % max voltage [p.u.]
Sinj = zeros(size(YBus,1),1); % inverter capacity vector 
Sinj(mpc.gen(2:end,1)) = sqrt(mpc.gen(2:end,9).^2+mpc.gen(2:end,4).^2); 
%%
Pav = zeros(size(YBus,1),1); % solar PV vector
Pav(mpc.gen(2:end,1)) = mpc.gen(2:end,9); % solar PV output [kW]
PF = 0.8; % selected power factor
theta = radtodeg(acos(PF)); 
Pd = Pd(2:end); % without slack bus
Qd = Qd(2:end); % without slack bus
%% Network visualisation
% network2 = digraph(linesFrom,linesTo, real(Z));  
% network2.Edges; % shows the number of Edges and Nodes 
% %A = adjacency(network); % shows all nodes that llnes are connected to 
% figure(2) 
% postNetwork.Edges.LWidths = 3*abs(network2.Edges.Weight)/max(abs(network2.Edges.Weight))
% if exist('mpc.bus_name','var') == 1
%     netgraph = plot(network2,'Layout','force','Linewidth',2,'NodeLabel',BusName,'ArrowSize',12,'EdgeLabel',round(network2.Edges.Weight,2));    
% else
%     netgraph = plot(network2,'Layout','force','Linewidth',2,'ArrowSize',12,'EdgeLabel',round(network2.Edges.Weight,2));        
% end
% highlight(netgraph,find(mpc.bus(:,2)==1),'NodeColor','black', 'MarkerSize',6); % loads 
% highlight(netgraph,find(mpc.bus(:,2)==3),'NodeColor','blue', 'MarkerSize',10); % gen 
% highlight(netgraph,find(mpc.bus(:,2)==2),'NodeColor','red', 'MarkerSize',6); % gen 
% set(gcf,'color','w'); 
%% Optimisation problem
cvx_begin 
    variable V(nBuses-1) complex; % Voltage vector [p.u.]
    variable Pc(nBuses-1); % Curtailed PV vector [kW]
    variable Pinj(nBuses-1); % Active power produced by inverter [kW]
    variable Qinj(nBuses-1); % Reative power absorbed/produced by inverter [kW] 
    Obj1 = 0;
    Objr = 0;
%     for m = 1 : nLines
%        for n = 1 : nLines
%            Obj1 = sum(real(YBus(m,n))*(real(V(m))+real(V(n)))^2); %Eq.1
%        end
%        Objr = Objr + Obj1;
%     end
    Obj2 = 0; 
    for n = 1 : size(YBus,1)
        % P+Q
        Obj2 = Obj2 + (250*Pc(n)+100*abs(Qinj(n))); % Eq.2 % 0.25$/kWh in MWh
        % P
        % Obj2 = Obj2 + (250*Pc(n)); % Eq.2 % 0.25$/kWh in MWh
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
       for n = 1 : nBuses-1
            Pav(n) - Pinj(n) <= Pc(n); % Eq.9
            0 <= Pinj(n) <= Pav(n); % Eq.10  
            square(Qinj(n)+Pinj(n)) <= (Sinj(n))^2; % Eq.11
            abs(Qinj(n)) <= tan(theta*pi()/180)*Pinj(n); % Eq.12
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
       % Min & Max voltage magnitude limit 
           for l = 1 : size(YBus,1)
               sub_minV = sub_minV + real(ZBus(n,l))*(Pinj(l) - Pd(l)) + imag(ZBus(n,l))*(Qinj(l) - Qd(l));
               sub_maxV = sub_maxV + real(ZBus(n,l))*(Pinj(l) - Pd(l)) + imag(ZBus(n,l))*(Qinj(l) - Qd(l)); 
           end
           Vmin(n) <= Vnom(n) + sub_minV; % Eq.7 
           Vnom(n) + sub_maxV <= Vmax(n); % Eq.8
       % Max voltage magnitude limit 
           %Vnom(n) + sum(real(ZBus(n,:))*(Pinj(n) - Pd(n))) + sum(imag(ZBus(n,:))*(Qinj(n) - Qd(n))) <= Vmax(n); % Eq.8  
       end  
cvx_end
% Results
Guggilam_V = [Vnom(1); real(V)]; % adding slack bus
if exist('mpc.bus_name','var') == 1
    VoltageTable = table(Guggilam_V, ACOPF_V, 'rownames', BusName) % compare Guggilam and ACOPF
    VoltageTable.Properties.VariableNames = {'Guggilam' 'ACOPF'};
    writetable(VoltageTable, 'V_Guggilam_VS_ACOPF', 'WriteRowNames',true);
    %
    ParTable = table(Pav, round(Pc,4), round(Pinj,4), round(Qinj,4), round(Sinj,4), 'rownames', BusName(2:end)); % compare Guggilam and ACOPF
    ParTable.Properties.VariableNames = {'Pav', 'Pc', 'Pinj', 'Qinj','Sinj'}
    writetable(ParTable, 'Par_Guggilam', 'WriteRowNames',true);
else 
    VoltageTable = table(Guggilam_V, ACOPF_V) % compare Guggilam and ACOPF
    VoltageTable.Properties.VariableNames = {'Guggilam' 'ACOPF'};
    writetable(VoltageTable, 'V_Guggilam_VS_ACOPF', 'WriteRowNames',false);
    %
    ParTable = table(Pav, round(Pc,4), round(Pinj,4), round(Qinj,4), round(Sinj,4)); % compare Guggilam and ACOPF
    ParTable.Properties.VariableNames = {'Pav', 'Pc', 'Pinj', 'Qinj','Sinj'}
    writetable(ParTable, 'Par_Guggilam', 'WriteRowNames',false);
end
% plot voltages 
figure(1)
plot(1:nBuses,[real(Vnom(1)); real(V)])
hold on
plot(1:nBuses,ACOPF_V)
xlabel('bus') 
ylabel('Voltage [p.u.]') 
xlim([1 nBuses]); ylim([0.9 1.1]) 
xticks(1:1:nBuses); 
legend('Guggilam', 'ACOPF')
set(gcf,'color','w'); 
