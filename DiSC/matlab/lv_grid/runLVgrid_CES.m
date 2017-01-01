 % Script for simulating the LV benchmark grid
clc;
clear all; close all;


% Call setup script
LVgrid_setup

% Allocate memory
Vout = zeros(numSamples,numBus);
Pslack = zeros(numSamples,1);
Qslack = zeros(numSamples,1);

% Control variables
bufferLength = 10;
maxDelay = 5;
Vbuf = ones(bufferLength,numBus);
Tbuf = zeros(bufferLength,numBus);
Qrefbuf = ones(bufferLength,numBus);
Trefbuf = zeros(bufferLength,numBus);

% Run power flow
for i=1:numSamples
    out = nrLoadFlow(Y,type,Pin(i,:)',Qin(i,:)',Vin',Vin',tol,maxIte);
    Vout(i,:) = out.Vout;
    Pslack(i) = out.Pslack;
    Qslack(i) = out.Qslack;
    if out.nrIte == 100
        disp('ERROR nrIte = 100')
        break
    end
    % Measurements to LVGC delay
    Vbuf(mod(i-1,bufferLength)+1,:) = Vout(i,:);
    Tbuf(mod(i-1,bufferLength)+1,:) = i+ceil(maxDelay*rand(1,numBus));
    
    % Current measurement
    [~,IndexMes] = max(Tbuf.*(Tbuf<=i),[],1);
    for j = numBus:-1:1
        Vmes(j) = Vbuf(IndexMes(j),j);
    end
    
    % LVGC
    alpha = zeros(numBus,numBus);
    alpha(32,32) = -1.5;
    alpha = -0.15*eye(numBus);
    Qsent = alpha*(abs(Vmes')-1);
    % Saturation
    
    % LVGC to household delay
    Qrefbuf(mod(i-1,bufferLength)+1,:) = Qsent;
    Trefbuf(mod(i-1,bufferLength)+1,:) = i+ceil(maxDelay*rand(1,numBus));
    
    % Current reference
    [~,IndexRef] = max(Trefbuf.*(Trefbuf<=i),[],1);
    for j = numBus:-1:1
        Qref(j) = Qrefbuf(IndexRef(j),j);
    end
    
    % Household
    if i<numSamples
        Qin(i+1,:) = Qin(i+1,:)+Qref;
    % Saturation
    end
    
end

%% Plotting
tvec = 0:numSamples-1;

figure
plot(tvec,Pslack*Sbase/1000,tvec,Qslack*Sbase/1000)
grid
title('Active and Reactive Consumption of Aggregated Grid')
legend('P','Q')
xlabel('Time [15 min]')
ylabel('Power [kW/kVAR]')

figure
plot(tvec,abs(Vout(:,15)),tvec,abs(Vout(:,32)),tvec,abs(Vout(:,42)))
grid
title('Voltage Profiles')
legend('bus15','bus32','bus42')
xlabel('Time [15 min]')
ylabel('Voltage p.u.')
