clc; close all;clear;

Ts = 10;
N = 200;
% System 
numP = 2;
A = [-1/(30*Ts) 0; 0 -1/(35*Ts)];
B = -A;
C = [1 0;0 1];
Cz = [1 1];
D = 0;

% Disturbance model
Ad = [0 0.01 0 0;-0.01 0 0 0;0 0 0 0.05; 0 0 -0.05 0];
Cd = [10 0 20 0];

% loop with disy
Al = [A B*((ones(numP,1)/numP)*Cd); zeros(4,2) Ad];
Bl = [B;zeros(4,2)];
Cl = [C zeros(2,4)];
sysl = ss(Al,Bl,Cl,0);
sysld = c2d(sysl,Ts);

sys = ss(A,B,C,D);
sysd = c2d(sys,Ts);

Q = eye(2)*1;
R = eye(2)*1;

Kd = dlqr(sysd.a,sysd.b,Q,R);
Ld = dlqr(sysld.a',sysld.c',eye(6),eye(2));
Ld = Ld'

% Augmented system
AI = [A zeros(2,1);Cz 0];
BI = [B;0 0];
CI = [Cz 0];
DI = 0;

sysI = ss(AI,BI,CI,DI);
sysdI = c2d(sysI,Ts)
KI = dlqr(sysdI.a,sysdI.b,eye(3),eye(2))


x = zeros(2,N);
y = zeros(2,N);
u = zeros(2,N);
x(:,1) = [10;5];
xd = zeros(4,N);
xd(:,1) = 1;
yd = zeros(1,N);

xhat = zeros(2,N);
yhat = zeros(2,N);
e = zeros(1,N);
ref = 0;
for i = 1:N
    if i == 50
        ref = 100;
    end
    x(:,i+1) = sysd.a*x(:,i) + sysd.b*u(:,i) + sysd.b*(ones(numP,1)/numP)*yd(i);
    y(:,i) = sysd.c*x(:,i);
    xd(:,i+1) = Ad*xd(:,i);
    yd(i) = Cd*xd(:,i);
    
    xhat(:,i+1) = sysd.a*xhat(:,i) + sysd.b*u(:,i) + Ld(1:2,1:2)*(y(:,i)-sysd.c*xhat(:,i));
    e(i+1) = e(i) + (sum(y(:,i))-ref);
    
    u(:,i+1) = -KI(1:2,1:2)*xhat(:,i+1);% -KI(1:2,3)*e(i+1);
end

%% Plotting
t = 0:N;

figure
plot(t,x,t,xhat)

figure
plot(u')

figure
plot(sum(y))