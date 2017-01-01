% System of systems course example grid
clc; close all; clear;

% General setup
numSamples = 100;

% Line impedances
% Zij = Rij + I*Xij
R12 = 0.1; X12 = 0.1;
R23 = 0.1; X23 = 0.1;
R34 = 0.1; X34 = 0.1;
R45 = 0.1; X45 = 0.1;
R56 = 0.1; X56 = 0.1;
R67 = 0.1; X67 = 0.1;

% Construct grid
% Format: [From To R X l], l is length of cable
Z = [1 2 R12 X12 1;
     2 3 R23 X23 1;
     3 4 R34 X34 1;
     4 5 R45 X45 1;
     5 6 R56 X56 1;
     6 7 R67 X67 1;];
nBus = max(Z(:,2));
% Setup load flow object
param.type = [0 ones(1,nBus-1)];
param.vBase = 400*ones(1,nBus);

% Construct power flow object
pf = powerFlow(param);
 
% Get admittance matrix
Y = pf.lDataToY(Z);

% Setup busses
PF = 0.95;          % Power factor
Pbus1 = 0;
Pbus2 = -8e3;
Pbus3 = -8e3;
Pbus4 = -8e3;
Pbus5 = -8e3;
Pbus6 = -8e3;
Pbus7 = -8e3;

Pin = [Pbus1 Pbus2 Pbus3 Pbus4 Pbus5 Pbus6 Pbus7];

Qin=Pin.*tan(acos(PF));

Vin = [400 zeros(1,nBus-1)]; % Bus1 = slack

% Allocate memory
vOut = 400*ones(nBus,numSamples+1);
pSlack = zeros(1,numSamples);
qSlack = zeros(1,numSamples);

% Run simulation
for i = 1:numSamples
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Make changes to Pin and Qin by e.g. Control
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve load flow
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [vOut(:,i+1),pSlack(i),qSlack(i),~]=pf.nrLoadFlow(Y,Pin',Qin');
end

%% Plotting
tvec = 0:numSamples;

figure
plot(tvec,abs(vOut))
xlabel('Time [Sample]')
ylabel('Voltage Magnitude [V]')
grid