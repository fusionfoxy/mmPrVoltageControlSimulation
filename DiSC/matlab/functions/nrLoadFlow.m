function [vOut,pSlack,qSlack,qOut,ite] = nrLoadFlow(ybus,type,Pin,Qin,vBase,tol,maxIte)
% nrLoadFlow is a Newton-Raphson implementation of electrical grid Load Flow
% calculation.
%
%   Calculates complex bus voltages of a connected Power system, based on
%   the Newton-Raphson (N-R) method for solving the nonlinear algebraic
%   equations. The N-R method is based on polar coordinates. 
%
%   Input:
%       - ybus, is the bus admitance matrix of the system dim(n,n), where n
%         is number of busses.
%       - type, defines the type of each bus: 0=slack, 1=PQ, 2=PV. The
%         algorithm uses the type for easyer calculation. dim(n)
%       - Pin, is the real power at each bus. positive for
%         production and negative for consumption. dim(n)
%       - Qin, is the reactive power at each bus. positive for production and
%         negative for consumption. If V is specified at the bus then Q is
%         typically not needed. dim(n)
%       - Vin, is the voltage magnitude at PV busses. That is, where some unit is
%         keeping a constant voltage level, e.g., a CHP. dim(n)
%       - tol, is the tolerance of the algorithm typically small e.g. 0.0005.
%       - maxIte, is the maximum number of iterations.
%
%   Output:
%       - Vout, is the complex voltage level at each bus
%       - Pslack, is the real power at the slack busses.
%       - Qslack, is the reactive power at the slack busses.
%       - Qout, is the reactive power at each PV bus.
%
%
%   R. Pedersen 5-26-2014, Aalborg University

%% Error handling
% Check for wrong input type
if ~isempty(find(type>2, 1)) || ~isempty(find(type<0, 1))
    error('Wrong type input for one or more of the bus types')
end

%% Initialization
% Error and iterations
err = 1;
ite = 0;

% Setup index parameters
nBus = length(ybus(:,1));               % Number of busses
nSlack = nBus - length(find(type));     % Number of slack busses
nPV = length(find(type>1));             % Number of PV busses
nPQ = nBus-nSlack-nPV;                  % Number of PQ busses

idSlack = find(type==0);                % Index of slack busses
idPQ = find(type==1);                   % Index of PQ busses
idPV = find(type==2);                   % Index of PV busses
idPQV = find(type>0);                   % Index of both PQ and PV busses


% Admitance matrix in polar coordinates
Ym = abs(ybus);                         % Magnitude of Ybus 
Ya = angle(ybus);                       % Angle of Ybus

% Initial guess on voltage magnitudes and angles                           (STEP 1)
Vm = zeros(nBus,1);                     % Vector of bus voltage magnitudes  
Va = zeros(nBus,1);                     % Vector of bus voltage angles

% Set slack bus voltages in Vm and Va vector
for k=1:nSlack
    Vm(idSlack(k)) = abs(vBase(idSlack(k)));
    Va(idSlack(k)) = angle(vBase(idSlack(k)));
end
% Set PV bus voltages in Vm and Va vector 
for k=1:nPV
    Vm(idPV(k)) = abs(vBase(idPV(k)));
    Va(idPV(k)) = angle(vBase(idSlack(1)));
end
% set PQ voltages in Vm and Va vector set to Vbase
for k=1:nPQ
    Vm(idPQ(k)) = abs(vBase(idPQ(k)));
    Va(idPQ(k)) = angle(vBase(idPQ(k)));
end

%% Run the Newton-Raphson Method
while err >= tol && ite < maxIte
    
    % Iterate
    ite = ite+1;
    
    % Allocate Memory
    P = zeros(nBus,1);
    Q = zeros(nBus,1);
    J1 = zeros(nBus-nSlack,nBus-nSlack);
    J2 = zeros(nBus-nSlack,nPQ);
    J3 = zeros(nPQ,nBus-nSlack);
    J4 = zeros(nPQ,nPQ);
    
    % For load busses (PQ), calculate P, Q.                                (STEP 2)
    for k=1:nPQ
        for n=1:nBus
            P(idPQ(k)) = P(idPQ(k)) + Vm(idPQ(k))*Vm(n)*Ym(idPQ(k),n) * ...
                            cos(Ya(idPQ(k),n)-Va(idPQ(k))+Va(n));
            Q(idPQ(k)) = Q(idPQ(k)) - (Vm(idPQ(k))*Vm(n)*Ym(idPQ(k),n) * ...
                            sin(Ya(idPQ(k),n)-Va(idPQ(k))+Va(n)));
        end
    end
    
    % For voltage-controled busses (PV), calculate P.                      (STEP 3)
    for k=1:nPV
        for n=1:nBus
            P(idPV(k)) = P(idPV(k)) + Vm(idPV(k))*Vm(n)*Ym(idPV(k),n) * ...
                            cos(Ya(idPV(k),n)-Va(idPV(k))+Va(n));
        end
    end
    
    % Calculate the elements J1,J2,J3 and J4 of the Jacobian matrix        (STEP 4)
    % J1 has dim(nBus-nSlack,nBus-nSlack)
    for k=1:nBus-nSlack
        % Diagonal elements of J1
        for n=1:nBus
            if idPQV(k)==n
            else
                J1(k,k) = J1(k,k) + Vm(idPQV(k))*Vm(n)*Ym(idPQV(k),n) * ...
                            sin(Ya(idPQV(k),n)-Va(idPQV(k))+Va(n));
                        
            end
        end
        % Off diagonal elements of J1
        for n=1:nBus-nSlack
           if k==n
           else
                J1(k,n) = -(Vm(idPQV(k))*Vm(idPQV(n))*Ym(idPQV(k),idPQV(n)) * ...
                            sin(Ya(idPQV(k),idPQV(n))-Va(idPQV(k))+Va(idPQV(n))));
           end
       end
    end
    
    % J2 has dim(nBus-nSlack,nPQ)
    for k=1:nPQ
        % Diagonal elements of J2
        for n=1:nBus
            if idPQV(k)==n
            else
                J2(k,k) = J2(k,k) + Vm(n)*Ym(idPQV(k),n) * ...
                             cos(Ya(idPQV(k),n)-Va(idPQV(k))+Va(n));
            end
        end
        % Add last element to diagonal of J2
        J2(k,k) = J2(k,k) + 2*Vm(idPQV(k))*Ym(idPQV(k),idPQV(k)) * ... 
                    cos(Ya(idPQV(k),idPQV(k)));
        % Off diagonal elements of J2
        for n=1:nBus-nSlack
            if k==n
            else
               J2(n,k) = Vm(idPQV(n))*Ym(idPQV(n),idPQV(k)) * ...
                            cos(Ya(idPQV(n),idPQV(k))-Va(idPQV(n))+Va(idPQV(k))); 
            end
        end
    end
    
    % J3 has dim(nPQ,nBus-nSlack)
    for k=1:nPQ
        % Diagonal elements of J3
        for n=1:nBus
            if idPQV(k)==n
            else
                J3(k,k) = J3(k,k) + Vm(idPQV(k))*Vm(n)*Ym(idPQV(k),n) * ...
                            cos(Ya(idPQV(k),n)-Va(idPQV(k))+Va(n));
            end 
        end
        % Off diagonal elements of J3
        for n=1:nBus-nSlack
            if k==n
            else
                J3(k,n) = -(Vm(idPQV(k))*Vm(idPQV(n))*Ym(idPQV(k),idPQV(n)) * ...
                            cos(Ya(idPQV(k),idPQV(n))-Va(idPQV(k))+Va(idPQV(n))));
            end
        end
    end
    
    % J4 has dim(nPQ,nPQ)
    for k=1:nPQ
        % Diagonal elements of J4
        for n=1:nBus
            if idPQV(k)==n
            else
                J4(k,k) = J4(k,k) - (Vm(n)*Ym(idPQV(k),n) * ...
                            sin(Ya(idPQV(k),n)-Va(idPQV(k))+Va(n)));
            end
        end
        % Add last element to diagonal of J2
        J4(k,k) = J4(k,k) - 2*Vm(idPQV(k))*Ym(idPQV(k),idPQV(k)) * ...
                    sin(Ya(idPQV(k),idPQV(k)));
        % Off diagonal elements of J4
        for n=1:nPQ
            if k==n
            else
                J4(k,n) = -(Vm(idPQV(k))*Ym(idPQV(k),idPQV(n)) * ...
                            sin(Ya(idPQV(k),idPQV(n))-Va(idPQV(k))+Va(idPQV(n))));
            end
        end
    end
    
    % Form J
    J = [J1 J2; J3 J4];
    
    % Find residuals and solve linear simultaneous equations               (STEP 5)
    dP = Pin(idPQV) - P(idPQV);
    dQ = Qin(idPQ) - Q(idPQ);
    dPQ = [dP; dQ];
    
    dVam = J\dPQ;
    
    % Compute new voltage magnitudes and phase angles                      (STEP 6)
    Va(idPQV) = Va(idPQV) + dVam(1:length(Va(idPQV)));
    Vm(idPQ) = Vm(idPQ) + dVam(length(Va(idPQV))+1:end);
    
    % Calculate error                                                      (STEP 7)
    err = max(abs(dPQ));
    
end

%% Set outputs
% Voltages at each bus including slack bus
vOut = Vm.*exp(1i*Va);

% Active and reactive power at slack bus
pSlack = zeros(nSlack,1);
qSlack = zeros(nSlack,1);
for k=1:nSlack
    for n=1:nBus
       pSlack(k) = pSlack(k) + Vm(k)*Vm(n)*Ym(k,n) * ...
                    cos(Ya(k,n)-Va(k)+Va(n));
       qSlack(k) = qSlack(k) - (Vm(k)*Vm(n)*Ym(k,n) * ...
                    sin(Ya(k,n)-Va(k)+Va(n))); 
    end
    
end

% Reactive power at PV busses
qOut = zeros(nPV,1);
for k=1:nPV
    for n=1:nBus
        qOut(k) = qOut(k) - (Vm(idPV(k))*Vm(n)*Ym(idPV(k),n) * ...
                    sin(Ya(idPV(k),n)-Va(idPV(k))+Va(n))); 
    end
end
    
end

