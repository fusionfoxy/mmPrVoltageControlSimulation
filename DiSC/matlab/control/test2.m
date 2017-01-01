clc; close all; clear;

Ts = 10;
N=ceil(1*24*60*60/Ts);
numLoad = 1;
numProd = 1;
numOsci = 2;
numXhat = numLoad*numOsci*2+1;
numYhat = 1;

% Disturbance model
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
A = [0      beta1   0       0;
     -beta1 0       0       0;
     0      0       0       beta2;
     0      0       -beta2  0;];
% C = [a1*cos(theta1) a1*sin(theta1) a2*cos(theta2) a2*sin(theta2) a3*cos(theta3) a3*sin(theta3)];   
C = [a1*cos(theta1) a1*sin(theta1) a2*cos(theta2) a2*sin(theta2)]; 
% B = [0;0;0;0;0;0];
B = [0;0;0;0];

Ad = blkdiag(A,0);
Bd = [B;0];
Cd = [C 1];
Dd = 0;

sys = ss(Ad,Bd,Cd,Dd);
sysdD = c2d(sys,Ts);

% System
% Controllable models
Ap = -1/(60);
Bp = -Ap;
Cp = 1;
Dp = 0;

AI = [Ap 0; Cp 0];
BI = [Bp;0];


sysp = ss(Ap,Bp,Cp,0);
syspd = c2d(sysp,Ts);

F = -lqr(Ap,Bp,eye(1),eye(1));
L = -lqr(Ap',Cp',eye(1),eye(1))';

% Zero assignment
Aza = Ap+Bp*F(1)+L*Cp;
Cza = -F(1);

Mt = -place(Aza',Cza',eig(Ap+Bp*F(1)))';

% Obtain unity DC-gain
Acl = [Ap Bp*F(1);-L*Cp Ap+Bp*F(1)+L*Cp];
Bclt = [
    Bp;Mt]; %(without N matrix)
Ccl = [Cp zeros(size(Cp,1),1)];

Nc = -eye(1)/(Ccl/Acl*Bclt);

M = Mt*Nc;

Bcl = Bclt*Nc;

SYS = ss(Acl,Bcl,Ccl,0);

figure
pzmap(SYS)

figure

bode(SYS)

figure
step(SYS)

% Controller
Ac = Ap+Bp*F(1)+L*Cp;
Bc = [M -L];
Cc = F(1);
Dc =[Nc 0];

Ctrl = ss(Ac,Bc,Cc,Dc);
% Accl = [A+B*Dc B*Cc;Bc*C Ac];
% Cccl = [C zeros(1,2);Dc*C Cc]

Ctrld = c2d(Ctrl,Ts);


%% Allocate memory
x = zeros(5,N);
ym = zeros(1,N);
ym2 = zeros(1,N);
ym3 = zeros(1,N);
x(:,1) = [1;0;1;0;b];
xhat = zeros(numLoad*numXhat,N);
xhat(:,1) = ones(numXhat,1);
yhat = zeros(numLoad*numYhat,N);
y = zeros(1,N);
u = zeros(1,N);
xP = zeros(1,N);
xP(1) = 500;
yP = zeros(1,N);
xPhat = zeros(1,N);
xc = zeros(1,N); 
e = zeros(1,N);


% For estimator
Q = eye(numXhat)*100;
Q(end,end) = 50;
P(:,:,1) = eye(numXhat);
R = 10000;
K = zeros(numXhat,N);

% Reference
ref = 0;

%% Simulation
for i=1:N
    % Disturbance
    x(:,i+1) = sysdD.a*x(:,i);
    y(i) = sysdD.c*x(:,i);
    ym(i) = y(i);
    
    % System
    xP(i+1) = syspd.a*xP(i) + syspd.b*u(i);
    yP(i) = syspd.c*xP(i);
    
    % Estimator (Disturbance)
    [xhat(:,i+1), yhat(i),phi,H] = periodicModel(numOsci,Ts,xhat(:,i));
    % Calculate covariance
    P(:,:,i+1) = phi*P(:,:,i)*phi' + Q;
    
    % Updating step
    % Calculate Kalman gain:
    K(:,i) = P(:,:,i+1)*H'/(H*P(:,:,i+1)*H'+R);
    
    % Update state estimate:
    xhat(:,i+1) = xhat(:,i+1)+K(:,i)*(ym(i)-yhat(i));
    
    % Update covariance estimate:
    P(:,:,i+1) = (eye(numXhat)-K(:,i)*H)*P(:,:,i+1);
    
    % Put in controller
    if i == 2000;
        ref = 100;
    elseif i > 3000;
        ref = -100;
    end
   
    e(i+1) = e(i) + ((ym(i)+yP(i))-ref)*Ts;
    xc(:,i+1) = Ctrld.a*xc(:,i) + Ctrld.b*[-(H*xhat(:,i+1)-ref);yP(i)];
    u(:,i+1) = Ctrld.c*xc(:,i+1) + Ctrld.d*[-(H*xhat(:,i+1)-ref);yP(i)] + (-0.01)*e(i);
end

%%
% xt = zeros(5,1440);
% yt = zeros(1,1440);
% xt(:,1) = xhat(:,end);
% for i=1:1440
%     xt(:,i+1) = phi*xt(:,i);
%     yt(i) = H*xt(:,i);
% end

%% Plotting
close all;
t=0:N-1;
figure
plot(t,ym,t,yhat,t,yP)
legend('ym','yhat','yP')

figure
plot(t,xhat(:,1:end-1))
legend('x1','x2','x3','x4','x5')

figure
plot(t,ym+yP)

% figure
% plot(t,K')