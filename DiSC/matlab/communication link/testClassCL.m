% Script for testing the communication link class
clc; close all; clear;

% Control variables
numLinks = 2;                           % Number of busses
N=10;                                   % Number of samples                      

% Setup communication link
param.Ts = 60;                          % Sampling time
param.maxDelay = 5*60;                  % Maximum delay
param.numLinks = numLinks;              % Number of communication links
param.mu = 2*60;
param.sigma = 1*60;
% Create communication link object
cl = comLink(param);
cl.setInverseCDFdelay('uniform',param);
cl.setPrLoss(0.3);

Vout = rand(N,numLinks);
%%
for i=1:N
    % Measurements to control communication
    Vmes = cl.sampleIn(i,Vout(i,:));
    
    % Controller
    alpha = -0.15*eye(numLinks);
    Qsent = alpha*(abs(Vmes')-1);
    % Saturation
    
    % Controller to asset communication
    Qref2=cl.sampleOut(i,Qsent);
end

%% Plotting
