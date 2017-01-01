clc
clear
close all

%
Ts = 60;
ite = 100;


A = [2 -3;4 -5];
B = [2 3]';
% B = diag(B);
C = [-3 2];
n = size(A,1);
p = size(C,1);
sysO = ss(A,B,C,0);

F = -place(A,B,[-3 -4]);
L = -place(A',C',[-9 -12])';

%% Zero assignment
Aza = A+B*F+L*C;
Cza = -F;

Mt = -place(Aza',Cza',[-3 -4])';

% Obtain unity DC-gain

Acl = [A B*F;-L*C A+B*F+L*C];
Bclt = [
    B;Mt]; %(without N matrix)
Ccl = [C zeros(size(C,1),n)];

% N = -eye(p)/(Ccl/Acl*Bclt);
N = -inv(Ccl/Acl*Bclt);

M = Mt*N;

Bcl = Bclt*N;

SYS = ss(Acl,Bcl,Ccl,0);

figure
pzmap(SYS)

figure

bode(SYS)

figure
step(SYS)




break
%% Collect control system
Ac = A+B*F+L*C;
Bc = [M -L];
Cc = F;
Dc =[N 0];

Ctrl = ss(Ac,Bc,Cc,Dc);
% Accl = [A+B*Dc B*Cc;Bc*C Ac];
% Cccl = [C zeros(1,2);Dc*C Cc]

Ctrld = c2d(Ctrl,Ts);

% Simulate
sys = ss(A,B,C,0);
sysd = c2d(sys,Ts);


x = zeros(n,ite);
x(:,1) = 0;
xc = zeros(n,ite);
y = zeros(1,ite);
u = zeros(1,ite);
e = zeros(1,ite);
ref = 1;
for i=1:ite
    % System
    x(:,i+1) = sysd.a*x(:,i) + sysd.b*u(i);
    y(i) = sysd.c*x(:,i+1);
    
    % Control
    e(i+1) = e(i) + (y(i)-ref)*Ts;
    xc(:,i+1) = Ctrld.a*xc(:,i) + Ctrld.b*[ref;y(i)];
    u(i+1) = Ctrld.c*xc(:,i+1) + Ctrld.d*[ref;y(i)];% + F(2)*e(i);
end

% Plotting
tvec = 0:ite-1;
figure
plot(tvec,y)