% Script for simulating the LV benchmark grid
% Needs: 
%       - LVgrid_setup.m
%       - LVlineImp.m
%       - nrLoarFlow.m
%       - ztoybus.m
%       - benchmarkImpedanceMatrix.m
%
% R. Pedersen 5-26-2014, Aalborg University
clc;
clear; close all;

Ts = 60;                   % Sampling time
% Call setup script
LVgrid_setup

% Allocate memory
Vout = ones(numSamples+1,numBus);
Pslack = zeros(numSamples,1);
Qslack = zeros(numSamples,1);

si = zeros(numSamples,numSI);
pPV = zeros(numSamples,numPV);
qPV = zeros(numSamples,numPV);
pES = zeros(numSamples,numES);
qES = zeros(numSamples,numES);
eES = zeros(numSamples,numES);

% Set input
day = 181;                  % start date: Julian date of the year (1-365)

dP = ones(1,numBus)*(0);     % [PU], active power change
dPlim = ones(1,numBus)*0.00;        % [PU], reference to derated power
qRef = -ones(1,numBus)*0.01;    % [PU], reactive power reference
vRef = ones(1,numBus);          % [PU], Voltage reference (for local droop ctrl.)
dP(esBus) = 0.01*0;

% Run simulation
tic
for i=1:numSamples
    % Itterate day
    if ~mod(i,24*60*60/Ts)
       if day == 365
           day = 1;
       else
           day = day +1;
       end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Insert production and flexible consumption
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if i==1
    else
        for n=1:numSI
            si(i,n) = siSys(n).sample(i,day,wEnvS(i));
        end
        for n=1:numPV
            [pPV(i,n), qPV(i,n)] = pvSys(n).sample(si(i,n),abs(Vout(i-1,pvBus(n))),dP(pvBus(n)),dPlim(pvBus(n)),qRef(pvBus(n)),vRef(pvBus(n)));
            Pin(i,pvBus(n)) = Pin(i,pvBus(n)) + pPV(i,n);
        end
        for n=1:numES
            [pES(i,n), qES(i,n), eES(i,n)] = esSys(n).sample(abs(Vout(i-1,esBus(n))),dP(esBus(n)),vRef(esBus(n)));
            Pin(i,esBus(n)) = Pin(i,esBus(n)) + pES(i,n);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Simulate Electrical Grid - LV  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    out=nrLoadFlow(Y,type,Pin(i,:)',Qin(i,:)',Vin',tol,maxIte);
    Vout(i,:) = out.Vout;       % Complex Voltage
    Pslack(i) = out.Pslack;
    Qslack(i) = out.Qslack;
    if out.nrIte == 100
        disp('ERROR nrIte = 100')
        break
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Measurement + Delay
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Apply Substation Control 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Send references + Delay
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Apply Local Control 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
end
toc
%% Plotting
tvec = (0:numSamples-1)/60;
tvecV = (0:numSamples)/60;

figure
plot(tvec,Pslack*Sbase/1000,tvec,Qslack*Sbase/1000)
grid
title('Active and Reactive Consumption of Aggregated Grid')
legend('P','Q')
xlabel('Time [hrs]')
ylabel('Power [kW/kVAR]')

figure
plot(tvecV,abs(Vout(:,15)),tvecV,abs(Vout(:,32)),tvecV,abs(Vout(:,42)))
grid
title('Voltage Profiles')
legend('bus15','bus32','bus42')
xlabel('Time [hrs]')
ylabel('Voltage p.u.')

figure
subplot(2,1,1)
plot(tvec,pPV*Sbase/1000)
grid
title('PV Systems')
ylabel('Active Power [kW]')

subplot(2,1,2)
plot(tvec,qPV*Sbase/1000)
grid
xlabel('Time [hrs]')
ylabel('Reactive Power [kVar]')

figure
subplot(3,1,1)
plot(tvec,pES*Sbase/1000)
grid
title('Energy Storages')
ylabel('Active Power [kW]')

subplot(3,1,2)
plot(tvec,qES*Sbase/1000)
grid
ylabel('Reactive Power [kVar]')

subplot(3,1,3)
plot(tvec,eES)
grid
xlabel('Time [hrs]')
ylabel('SOC [-]')