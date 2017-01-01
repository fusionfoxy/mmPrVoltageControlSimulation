clc
clear
close all

%% 

A = [2 -3;4 -5];
B = [2 3]';
C = [-3 2];
n = size(A,1);
p = size(C,1);

F = -place(A,B,[-3 -4]);
L = -place(A',C',[-9 -12])';

% Zero assignment
Aza = A+B*F+L*C;
Cza = -F;

Mt = -place(Aza',Cza',[-3 -4])';

% Obtain unity DC-gain

Acl = [A B*F;-L*C A+B*F+L*C];
Bclt = [
    B;Mt]; %(without N matrix)
Ccl = [C zeros(size(C,1),n)];

N = -eye(p)/(Ccl/Acl*Bclt);

%M = Mt*N;

Bcl = Bclt*N;

SYS = ss(Acl,Bcl,Ccl,0);

figure
pzmap(SYS)

figure

bode(SYS)

figure
step(SYS)
break
%%

Mt = -place(Aza',Cza',[-9 -12])';
Acl = [A B*F;-L*C A+B*F+L*C];
Bclt = [B;Mt]; %(without N matrix)
Ccl = [C zeros(size(C,1),n)];

N = -eye(p)/(Ccl/Acl*Bclt);
%M = Mt*N;
Bcl = Bclt*N;

SYS1 = ss(Acl,Bcl,Ccl,0);

%figure
hold on
%bode(SYS1)

% figure
% pzmap(SYS1)

figure
step(SYS1)

break
%%

N = 0;
M = L;

Bcl = [B*N;M];

SYS3 = ss(Acl,Bcl,Ccl,0);
hold on
bode(SYS3)
% [Y,T,X] = step(SYS3) 
% 
% figure
% plot(T,Y)
%% LQR

A = [2 -3;4 -5];
B = [2 3]';
C = [-3 2];
n = size(A,1);
p = size(C,1);

F1 = -lqr(A,B,800*(C')*C,1);

N1 = -eye(p)/(C/(A+B*F1)*B);
F2 = -place(A,B,[-4 -8]);
N2 = -eye(p)/(C/(A+B*F2)*B);

SYS1 = ss(A+B*F1,B*N1,C,0);
SYS2 = ss(A+B*F2,B*N2,C,0);
figure
step(SYS1)
hold on
step(SYS2)
xlim([0 1.5])
ylim([0 1.1])

