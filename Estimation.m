clear all; close all; clc;

filePrefix = '5sec_0985lim';%NewMM
delays = 1; %[0:0.05:0.2 0.3:0.1:0.8 1 2]; %[0:0.1:1 1.5:0.5:10]; %  [0:0.2:1 2:1:6]; %% 
seed = 1; %[1:2 4:6 31:32 34:36]; %1:60; %
% reactive = false;

dataRea = load(['results/' filePrefix '_delay' num2str(delays) '_reactive1_seed' num2str(seed) '_fullData.mat']);%res/loss/
dataPro = load(['results/' filePrefix '_delay' num2str(delays) '_reactive0_seed' num2str(seed) '_fullData.mat']);%res/loss/



for i = 2:length(dataPro.vOut)
    