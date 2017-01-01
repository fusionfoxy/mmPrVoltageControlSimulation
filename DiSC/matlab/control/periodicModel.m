function [x,y,A,C] = periodicModel(numOsci,Ts,xold)
%Function for simulating a system of harmonic oscillators
%   This function can be used for simulating an arbitrary number of
%   harmonic oscillators. The output is the sum of these systems.
%   Each oscillator is on the skew symetric form:
%   
%   dx = [0 beta; -beta 0]x
%   y = [a*cos(theta) a*sin(theta)] + b
%
%   where a is a gain, theta is the phase, b is a bias therm and beta represents the
%   frequency.
%
%   An arbetrary number of oscillators can be combined by placing each of
%   the above systems in the diagonal. This is further explained in the
%   note.
%
%   Input:
%       - numOsci, is the number of harmonic oscillators in the system.
%       - Ts, is the sample time, used for discritization.
%       - xold, is the statevector from last sample. OBS, when applying a
%       extended Kalman filter as intended, the state vector also containes
%       parameters of the system, such as a and theta.
%
%   Output:
%       - y, is the output of the system.
%       - A, is the Jacobian of the system function (System is linear, but needed for EKF)
%       - C, is the Jacobian of the output functions (-||-).
%       - x, is the new state vector.

% Frequencies. timeperiod must be a fraction of 24 hours
T = 24*60*60;
beta = zeros(numOsci,1);
beta(1) = 2*pi*1/T;
A = [0 beta(1); -beta(1) 0];
for i = 2:numOsci
    beta(i) = 2*pi*i/T;
    At = [0 beta(i); -beta(i) 0];
    A = blkdiag(A,At);
end
% Insert zeros for constant states to be estimated
% At = zeros(2,2);
% for i = 1:numOsci
%     A = blkdiag(A,At);
% end
% Finalyse system matrix by inserting for bias term
A = blkdiag(A,0);
A = expm(A*Ts); % Discritise
% Sample system
x = A*xold;

% Form output matrix
% C = zeros(1,numOsci*2+1);
% for i=1:numOsci*2
%     C(i) = x(numOsci*2+i);
% end
% C(end) = 1;
C = [ones(1,numOsci*2) 1];

y = C*xold;
end

