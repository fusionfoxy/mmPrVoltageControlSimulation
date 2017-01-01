% MV-Grid Benchmark
clc; clear all; %close all;

% Script for setting up grid
MVgrid_setup

% Run power flow for 24 hours
Vout = zeros(12,24);
Pslack = zeros(1,24);
Qslack = zeros(1,24);
for i=1:24
    out=nrLoadFlow(Y,type,Pin(i,:)',Qin(i,:)',Vin',tol,maxIte);
    Vout(:,i) = out.Vout;
    out.nrIte;
    Pslack(i) = out.Pslack;
    Qslack(i) = out.Qslack;
    if i==10 || i == 15
        Vin(1) = Vin(1)*1.0;
    end
end

%% Plotting
figure()
plot(0:23,Vout(5:11,:))
legend('Bus5','Bus6','Bus7','Bus8','Bus9','Bus10','Bus11')
grid
figure
plot(0:23,Pslack,0:23,Qslack)
legend('Ps','Qs')
grid
