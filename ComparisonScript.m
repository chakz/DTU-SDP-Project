% Comparison between Mosek, cvx and MatPower
close all; clear all; clc 
addpath(genpath('matpower6.0'))
addpath(genpath('mosek'))
mpc = IEEE_9BUS_Radial_modified;
%% AC OPF
mpopt = mpoption('model','AC', 'pf.tol', 1e-4,'opf.ac.solver','DEFAULT');
ACOPF = runopf(mpc,mpopt);
% Storing Voltage Magnitude values from AC OPF
ACOPF_V = ACOPF.bus(:,8); % excluding slack bus (saving results for comparison) 

%% Defining Network Topology
% We get generation data from the Generators script provided by MatPower
[genMatrix,nGen, genLoc, baseMVA, PMin, PMax, QMin, QMax, nBuses, busLoc, Vmin, Vmax, Pd, Qd] = generators(mpc);
% The Network data pertaining to impedances is provided from the Lines
% script
[LDCincidenceMat, linesMatFrom, linesMatTo, nLines, linesFrom, linesTo, R, X, B, Z, lineMaxFlow, OriginBusLoc] = lines(mpc);
% adds bus names if given in MatPower file 
if ~exist('mpc.bus_name','var')
    BusName = mpc.bus_name; 
    %BusInfo = [num2cell(busLoc'), num2cell(OriginBusLoc), BusName];
end

%% Bus Admittance Matrix
% The admittance matrix inclusive of the slack bus is provided by the
% script makeYbus
YBus_nplus1 = makeYbus(mpc); 
% We now want to remove the slack bus influence for this we will remove the
% slack bus from the column data and the row data
YBus = YBus_nplus1(2:end,2:end);
% The impedance matrix inclusive of the slack bus
ZBus_nplus1 = LDCincidenceMat*diag(Z);
ZBus = ZBus_nplus1(2:end,:);

%% Define parameters 
% nominal voltage vector [p.u.] This vector is a N*1 vector of the nominal
% voltage points
Vnom = ones(size(YBus,1),1); 
% The lower bounds on the voltage are provided as follows
Vmin = ones(size(YBus,1),1)*0.95; % min voltage [p.u.]
% The upper bound on the voltages are provied as follows
Vmax = ones(size(YBus,1),1)*1.05; % max voltage [p.u.]
% The net complex power injected into the bus 
Sinj = ones(size(YBus,1),1)*0.1; % inverter capacity [kW]
% The available power from a Solar panel
Pav = ones(size(YBus,1),1)*0; % solar PV output [kW]
% The power factor given the ratio between the real power / total power
PF = 0.8; % power factor 
% The active power demand specified in the matpower case file
Pd = Pd(2:end);
Qd = Qd(2:end);

%% Defining the Optimization Problem using the CVX solver
% Variable declaration 
% variables in cvx are defined as <variable name> (space) <size of
% variable>
cvx_begin 
    variable V(nBuses-1) complex; % Voltage vector [p.u.]
    variable Pc(nBuses-1); % Curtailed PV vector [kW]
    variable Qc(nBuses-1); % Curtailed PV reactive power vector [kW]
    variable Pinj(nBuses-1); % Active power produced by inverter [kW]
    variable Qinj(nBuses-1); % Reative power absorbed/produced by inverter [kW] 
    
    %Defining objective function for minimizing line losses % Equation 1
    Obj1 = 0;
     sub_Obj1_real = 0;
     sub_Obj1_imag = 0;
     store_Obj1_real = 0;
     for m = 1:nLines
          for n = 1:nLines
              real_voltage = [real(V(m)) real(V(n))];
              sum_real_voltage = sum(real_voltage,2);
              imag_voltage = [imag(V(m)) imag(V(n))];
              sum_imag_voltage = sum(imag_voltage,2);
              test_R = real(conj(YBus(m,n)));
              value_test_R = nonzeros(test_R);
              mag_value_test_R = double(value_test_R);
              sub_Obj1_real = (sum_real_voltage)' * mag_value_test_R * (sum_real_voltage);
              sub_Obj2_imag = (sum_imag_voltage)' * mag_value_test_R * (sum_imag_voltage);
          end
     end
     Obj1 = sub_Obj1_real + sub_Obj2_imag;
    %Obj1 = 0;
    %Obj1 = (sum(real(V(1))))'*real(conj(YBus(1,1)))*sum(real(V(1)));
    %Obj1 = sub_Obj1_real + sub_Obj1_imag;
%      for m = 1:nLines
%         for n = 1:nLines
%             real_voltages = [real(V(m)) real(V(n))];
%             sum_real_voltages = sum(real_voltages);
%             imag_voltages = [imag(V(m)) imag(V(n))];
%             sum_imag_voltages = sum(imag_voltages);
%             sub_obj1 = sum_real_voltages' * 0.8 * sum_real_voltages;
%             sub_obj2 = sum_imag_voltages' * 0.8 * sum_imag_voltages;
%             %sub_obj1 = quad_form(sum_real_voltages, 10);
%             %sub_obj1 = quad_form(real(V(m)) + real(V(n)), 10);
%             %square_sum_real_voltages = square(sum_real_voltages);
%             %square_sum_imag_voltages = square(sum_imag_voltages);
%             %square_total_voltages = [square_sum_real_voltages square_sum_imag_voltages];
%             %sum_square_total_voltages = sum(square_total_voltages);
%             %real_imp = nonzeros(real(conj(YBus(m,n))));
%             %Obj1 = Obj1 + real(conj(YBus(m,n))) * sum_square(real_voltages);
%             Obj1 = Obj1 + sub_obj1 + sub_obj2;
%         end
%     end
    
    
%     for m = 1:nLines
%         for n=1:nLines
%             square_real_voltage = (real(V(m)) + real(V(n)))^2;
%             square_real_voltage = real(V(m))^2 + real(V(n))^2 + 2*real(V(m))*real(V(n));                
%             square_imag_voltage = (imag(V(m)) + imag(V(n)))^2;
%             square_imag_voltage = imag(V(m))^2 + imag(V(n))^2 + 2*imag(V(m))*imag(V(n));  
%             Obj1 = Obj1 + real(conj(YBus(m,n))) * (square_real_voltage + square_imag_voltage);
%         end
%     end
%     for m = 1 : nLines
%       for n = 1 : nLines
%           %square_real_voltage = real(V(m))^2 + real(V(n))^2 + 2*real(V(m))*real(V(n));
%           %square_imag_voltage = imag(V(m))^2 + imag(V(n))^2 + 2*imag(V(m))*imag(V(n));
%           sum_real_voltage = real(V(m)) + real(V(n));
%           sum_imag_voltage = imag(V(m)) + imag(V(n));
%           Obj1 = Obj1 + real(conj(YBus(m,n)))* (square(sum_real_voltage) + square(sum_imag_voltage)); %Eq.1
%       end
%     end

    
    % Defining the Objective for minimizing active power curtailment %
    % Equation 2
    Obj2 = 0;
    % Defining the coefficients associated with this objective function
     a_h = 0; % Quadratic cost term (active power)
     b_h = 80; % Linear cost term (active power)
     c_h = 0; % Quadratic cost term (reactive power)
     d_h = 0; % Linear cost term
    for n = 1:size(YBus,1)
        Obj2 = Obj2 + a_h * Pc(n)^2 + (b_h * Pc(n)); %+ c_h * Qc(n)^2 + d_h*abs(Qc(n));       
    end
    
    % Defining the objective function pertaining to voltage regulation %
    % Equation 3
    
    Obj3 = 0;
%     for n = 1:size(YBus, 2)
%         for l = 1:size(YBus,1)
%             neighbor_nodes = sum(abs(V(l)));
%             %neigbors = [abs(V(n)) neighbor_nodes];
%             %ref = sum(neigbors);
%             %Obj3 = Obj3 + (abs(V(n)) - (1/((nBuses-1)+1)*sum(abs(V(l)))))^2;
%             Obj3 = Obj3 + neighbor_nodes;
%         end 
%     end
    %temp = sum(abs(V*0.1));
    %Obj3 = sum(abs(V(:))-temp);
    %Obj3 = temp;
    % Adding everything together and expressing the objective function 
    minimize (Obj1 + Obj2 + Obj3) % Equation 4
    subject to
    % Adding the PV Constraints
    for n=1:size(YBus,2)
        Pav(n) - Pinj(n) == Pc(n); % Equation 9 which is an extra constraint
        0 <= Pinj(n) <= Pav(n); % Equation 10 these constraints do not make complete sense
        Qinj(n)^2 + Pinj(n)^2 <= Sinj(n)^2; %Equation 11
        Qinj(n) <= PF * Pinj(n); % Equation 12
    end
    
    % Adding constraints related to voltage values
    % This is for storing results of voltage magnitudes for each node
    store_voltage_real = 0;
    store_voltage_imag = 0;
    % This is for capturing the min and max values
    store_voltage_min = 0;
    store_voltage_max = 0;
    
    for n = 1:size(YBus,2)
        for l = 1:size(YBus,1)
            store_voltage_real = store_voltage_real + (real(ZBus(n,l))*(Pinj(l) - Pd(l))) + ...
                                 (imag(ZBus(n,l)) * (Qinj(l) - Qd(l)));
                             
            store_voltage_imag = store_voltage_imag + (imag(ZBus(n,l)) * (Pinj(l) - Pd(l))) - ...
                                 (real(ZBus(n,l)) * (Qinj(l) - Qd(l)));
        end
        real(V(n)) == Vnom(n) + store_voltage_real; % Equation 4
        % Not sure if this value should be the following
        % real(V(n)) == abs(Vnom(n)) + store_voltage_real;
        imag(V(n)) == Vnom(n) + store_voltage_imag; % Equation 5
        % This could be
        % imag(V(n)) == angle(Vnom(n)) + store_voltage_imag
        
        for l = 1:size(YBus, 1)
            store_voltage_min = store_voltage_min + (real(ZBus(n,l))*(Pinj(l) - Pd(l))) + ...
                                (imag(ZBus(n,l)) * (Qinj(l) - Qd(l)));
            store_voltage_max = store_voltage_max + (real(ZBus(n,l))*(Pinj(l) - Pd(l))) + ...
                                (imag(ZBus(n,l))*(Qinj(l) - Qd(l)));
        end
        Vmin(n) <= Vnom(n) + store_voltage_min;
        Vnom(n) + store_voltage_max <= Vmax(n);
    end
cvx_end

% Begin MOSEK code




% Comparing the results
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
    


    
    
    
    
    