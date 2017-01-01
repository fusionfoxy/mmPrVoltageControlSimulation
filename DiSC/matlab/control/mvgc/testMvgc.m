% Script for testing the MVGC
clc; close all; clear;
Ts = 10;
N=2*24*60*60/(Ts);

% Set parameters
numLoad = 3;
numProd = 2;
param.numLoad = numLoad;
%param.aProd = [-1/(30*Ts) -1/(35*Ts)];
param.aProd = [-1/(10) -1/(5)];
param.bProd = -param.aProd;
param.cProd = [1 1];
param.Ts = Ts;

% Create object
mvgc = mvgcM1(param);
% Tune observer
Qd = eye(param.numLoad*(2*4+1))*10;
% Qd(2*4+1:2*4+1:end,2*4+1:2*4+1:end) = 0.1;
% Qd(2*4+1:2*4+1:end,2*4+1:2*4+1:end) = 0.1;
Rd = eye(1)*10;
mvgc.tuneDistEKF(Qd,Rd);
% Tune controller
Qs = eye(3);
Qs(3,3) = 1;
Rs = eye(2);
%mvgc.tuneCtrl(Qs,Rs)


beta1 = 7.2694e-5;
beta2 = 1.4544e-4;
beta3 = 2.1819e-4;
a1 = 31.89430610;
a2 = 30.77542974;
a3 = 9.088672378;
theta1 = 1.624103808;
theta2 = 2.350402122;
theta3 = 1.227124039;
b = 97.3411;
%%
% Disturbance model
% A = [0      beta1   0       0       0       0;
%      -beta1 0       0       0       0       0;
%      0      0       0       beta2   0       0;
%      0      0       -beta2  0       0       0;
%      0      0       0       0       0       beta3;
%      0      0       0       0       -beta3  0];
A = [0      beta1   0       0;
     -beta1 0       0       0;
     0      0       0       beta2;
     0      0       -beta2  0;];
% C = [a1*cos(theta1) a1*sin(theta1) a2*cos(theta2) a2*sin(theta2) a3*cos(theta3) a3*sin(theta3)];   
C = [a1*cos(theta1) a1*sin(theta1) a2*cos(theta2) a2*sin(theta2)]; 
% B = [0;0;0;0;0;0];
B = [0;0;0;0];
D = 0;

sys = ss(A,B,C,D);
sysdD = c2d(sys,Ts);

% Controllable models
Ac = diag(param.aProd);
Bc = diag(param.bProd);
Cc = diag(param.cProd);
Dc = 0;

sysP = ss(Ac,Bc,Cc,Dc);
sysPd = c2d(sysP,Ts);

% Allocate memory
x = zeros(4,N);
ym1 = zeros(1,N);
ym2 = zeros(1,N);
ym3 = zeros(1,N);
x(:,1) = [1;0;1;0];
xhat = zeros(numLoad*(2*4+1),N);
yhat = zeros(numLoad,N);
y = zeros(1,N);
u = zeros(2,N);
xP = zeros(2,N);
xP(:,1) = [10;5];
yP = zeros(2,N);
yAll = zeros(numProd+numLoad,N);

ref = 0;
%%
for i = 1:N
    % Disturbance
    x(:,i+1) = sysdD.a*x(:,i);
    y(i) = sysdD.c*x(:,i) + b;
    ym1(i) = y(i);
    ym2(i) = y(i) - 20;
    ym3(i) = y(i) - 50;
    
    % Controllable units
    xP(:,i+1) = sysPd.a*xP(:,i) + sysPd.b*u(:,i);
    yP(:,i) = sysPd.c*xP(:,i);
    yAll(:,i) = [yP(:,i);ym1(i);ym2(i);ym3(i)];
    
    
    if i == 3000
        ref = 100;
    end
    [yhat(:,i),xhat(:,i),u(:,i+1)]= mvgc.sample(ref,yAll(:,i));
end

%% Plotting
close all;
t = 0:N-1;
figure
plot(t,ym1,t,ym2,t,ym3,t,yhat)
legend('ym1','ym2','ym3','yhat1','yhat2','yhat3')

figure
plot(t,xhat(:,:))

figure
plot(t,[ym1;ym2;ym3],t,sum([ym1;ym2;ym3]),t,sum(yP))
legend('ym1','ym2','ym3','sum(ym)','sum(yP)')

figure
plot(t,sum(yAll))