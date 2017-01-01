classdef mvgcM1 < handle
    %MVGCM1 (medium voltage grid controller - Mark 1) is a medium voltage
    % controller class for the simulation framework DiSC.
    %
    %   This class will implement a specific controller class for a medium
    %   voltage grid. Its assignments are to ensure voltages at each bus
    %   and make the entire MV grid follow an energy reference. Following
    %   an energy reference can be used for e.g., engaging in the energy
    %   markets or offering frequency control services. In practise the
    %   controller is designed to follow a power reference.
    %
    %   The controller also incorporates an observer for estimating the
    %   power consumption of uncontrollable loads and production from
    %   controllable production units. Controllable production units could
    %   be a wind or solar power plant. To estimate uncontrollable 
    %   aggregated consumption, oscillator models are used.
    %
    %   Production units are modelled by first order systems:
    %       dx(t)/dt = ax(t) + bu(t),
    %       y(t) = cx(t)
    %   where x is the state, a, b and c are parameters and u is the input,
    %   e.g., reference for active power. Depending on sample time the
    %   production units can be seen as static, i.e., y(t) = u(t).
    %
    %   Uncontrollable aggregated consumption is modelled by ocillators:
    %       dx(t/dt) = [0 beta; -beta 0]x(t)
    %       y(t) = [a*cos(theta) a*sin(theta)]x(t) + b
    %   where a is a gain, theta is the phase, b is a bias therm and beta
    %   represents the frequency. Each consumption model can consist of
    %   several oscillator models.
    
    % Parameters for the controller object
    properties
        % General parameters
        numProd = 0;        % [-]. Number of flexible units.
        numLoad = 0;        % [-]. Number of loads.
        numOsci = 2;        % [-]. Number of oscillators for each load model.
        Ts = 60;            % [s]. Sampling time.
        busType             % [-]. Type of each bus in the system (0 = free, 1 = load, 2 = Flexible).
        
        % Observer parameters
        Qd                  % [-]. State covariance of disturbance, can be used for tuning.
        Rd                  % [-]. Output covariance matrix, ca be used for tuning.
        PdOld               % [-]. Old state covariance, used for the Kalman filter.
        xhatOld             % [-]. Holds the old estimates of states.
        numDstates          % [-]. Number of disturbance states.
        numCstates          % [-]. Number of controllable states.
        
        % Control parameters
        uOld                % [-]. Previous value of the controller output.
        Ctrld               % [-]. Control state space system (discrete time).
        xcOld               % [-]. Previous control state.
        uIOld = 0;          % [-]. Control integral output.
        tauMin = 10;        % [-]. Minimum time constant for deciding to run with dynamic or static control law. If below tauMin it uses dynamic feedback
        wInt = false        % [-]. Indicates if the integral term of the controller is active.
        wDispatch = false   % [-]. Indicates if the controller should use optimal dispatch of active power.
        satFlagUp           % [-]. Flags indicating if assets are saturated upwards.
        satFlagDown         % [-]. Flags indicating if assets are saturated downwards.
        kI = -0.08          % [-]. Integral gain.
        kP = 1;             % [-]. Proportional gain.
        
        
        % Dispatch parameters
        Z                   % [Ohm]. Is the systems bus impedance matrix, used for loss minimization.
        numBus              % [-]. Number of busses.
        idFree              % [-]. Index of busses without any production or consumption.
        idLoad              % [-]. Index of load busses.
        idProd              % [-]. Index of production busses.
        PpOld               % [W]. Previous values of dispatched active power values.
    end
    
    methods
        % Constructor
        function obj = mvgcM1(param)
            % param format:
            % - param.aProd         [-]. Is a vector containing the 'a' parameter of each flexible unit
            % - param.bProd         [-]. Is a vector containing the 'b' parameter of each flexible unit
            % - param.cProd         [-]. Is a vector containing the 'c' parameter of each flexible unit
            % - param.numLoad       [-]. Number of uncontrollable aggregated loads
            % - param.Ts            [s]. Sampling time
            
            % Set general parameters
            obj.numProd = length(param.aProd);
            obj.numLoad = param.numLoad;
            obj.Z = param.Z;
            obj.Ts = param.Ts;
            obj.busType = param.busType;
            obj.numDstates = obj.numLoad*(obj.numOsci*4+1);
            obj.numCstates = obj.numProd;
            
            % Set parameters for observer
            % Disturbance estimator
            obj.Qd = eye(obj.numDstates);                   % State covariance (Tuning parameter)
            obj.Rd = eye(obj.numLoad);                      % Output covariance (Tuning parameter)
            obj.PdOld = zeros(obj.numDstates);              % Kalman filter covariance
            obj.xhatOld = ones(obj.numCstates+obj.numDstates,1);
            
            % Set parameters controller
            obj.uOld = zeros(obj.numProd,1);
            obj.uIOld = 0;
            obj.xcOld = zeros(obj.numProd,1);
            obj.tauMin = min(abs(1./param.aProd))/15;
            obj.satFlagUp = zeros(obj.numProd,1);
            obj.satFlagDown = zeros(obj.numProd,1);
            % Design controller without disturbance and observer for
            % controllable part.
            % Create system matrices
            As = diag(param.aProd);
%             Bs = diag(param.bProd);
            Bs = param.bProd';
            Cs = param.cProd;
            % Observer for controllable part
            L = -lqr(As',Cs',eye(length(As)),eye(size(Cs',2)))';
            % Calculate gain matrix based on lqr and sampling time
            K = -lqr(As,Bs,eye(length(As)),eye(size(Bs,2)));
            
            % Set parameters for feed-forward
            % Zero assignment
            Aza = As+Bs*K+L*Cs;
            Cza = -K;
            Mt = -place(Aza',Cza',eig(As+Bs*K))';
            
            Acl = [As Bs*K;-L*Cs As+Bs*K+L*Cs];
            Bclt = [Bs;Mt]; %(without N matrix)
%             Ccl = [Cs zeros(size(Cs,1),obj.numProd)];
            Ccl = [Cs zeros(size(Cs,1),size(As,1))];
            N = -eye(size(Cs,1))/(Ccl/Acl*Bclt);
            M = Mt*N;
            
            % Form control system
            Ac = As+Bs*K+L*Cs;
            Bc = [M -L];
            Cc = K;
            Dc =[N zeros(size(Cc,1))];
            Ctrl = ss(Ac,Bc,Cc,Dc);
            obj.Ctrld = c2d(Ctrl,obj.Ts);
            
            % Set parameters for dispatch
            obj.numBus = length(obj.busType);
            obj.idFree = find(obj.busType==0);
            obj.idLoad = find(obj.busType==1);
            obj.idProd = find(obj.busType==2);
            obj.PpOld = zeros(obj.numProd,1);
        end
        
        % Sample the controller
        function [yhat,xhat,uRef]= sample(obj,pRef,pAsset,pMeas,pMax,pMin,vBus)
                      
            % Estimation
            [xhat,yhat,zhat,Ad,Cd] = observer(obj,pAsset(obj.numProd+1:end));
            
             % Control
             uRef = controller(obj,pRef,pAsset(1:obj.numProd),pMeas,zhat,yhat,pMax,pMin,vBus);
        end
        
        % Observer
        function [xhat,yhat,zhat,Ad,Cd] = observer(obj,y)
            % The observer implements a Kalman filter for
            % estimating the states of the uncontrollable
            % disturbance. 
            %
            % Input:
            %   - y [-], is a vector of measurements of uncontrollable part.
            %
            % Output:
            %   - xhat [-], is the estimated states.
            %   - yhat [-], is the estimated output value of disturbance.
            %   - zhat [-], is the output to be controlled.
            %   - Ad [-], is the disturbance system matrix, updated as some parameters are estimated.
            %   - Cd [-], is the disturbance output matrix, updated as some parameters are estimated.
           
            % Sample periodic models
            Ad = zeros(obj.numDstates,obj.numDstates);
            Cd = zeros(obj.numLoad,obj.numDstates);
            Cdz = zeros(1,obj.numDstates);
            xhatd = zeros(obj.numDstates,1);
            yhatd = zeros(obj.numLoad,1);
            numOstates = obj.numOsci*4+1;
            for i=1:obj.numLoad
                [xhatd((i-1)*numOstates+1:i*numOstates), yhatd(i), A, C] = periodicModel(obj.numOsci,obj.Ts,obj.xhatOld((i-1)*numOstates+1:i*numOstates));
                Ad((i-1)*numOstates+1:i*numOstates,(i-1)*numOstates+1:i*numOstates) = A;
                Cd(i,(i-1)*numOstates+1:i*numOstates) = C;
                Cdz(1,(i-1)*numOstates+1:i*numOstates) = C;
            end
            % Calculate covariance
            P = Ad*obj.PdOld*Ad' + obj.Qd;
            % Calculate Kalman Gain
            K = P*Cd'/(Cd*P*Cd'+obj.Rd);
            % Update state estimate:
            xhatd = xhatd+K*(y-yhatd);
            % Update covariance estimate:
            P = (eye(obj.numDstates)-K*Cd)*P;

            % Collect output
            xhat = xhatd;
            yhat = Cd*xhatd;
            zhat = Cdz*xhatd;
            % Save states and covariance
            obj.xhatOld = xhat;
            obj.PdOld = P;
        end
         
        % Controller
        function u = controller(obj,ref,yProd,yMeas,zhat,yhat,pMax,pMin,vBus)
            % The controller is implemented as an internal model controller
            % which uses the disturbance model in the control apporach.
            % Further, the system has been augmented with an additional
            % integral state. The integral controller uses the measurements
            % from the controllable part directly, to not have any bias
            % which may occure if only using the state estimates. 
            %
            % Input:
            %   - ref [W], is the active power reference signal for the controller
            %   - yProd [W], is a vector of active power measurements at the production units
            %   - yMeas [W], is the measurement of the entire systems power (used for integral control)
            %   - zhat [W], is the disturbance estimate obtained by the observer
            %   - yhat [W], is a vector containing the estimated power consumption of each disturbance model
            %   - pMax [W], is a vector containing information on the maximum power available of the flexible assets (used for integral anti-windup and dispatch)
            %   - pMin [W], is a vector containing information on the minimum power available of the flexible assets (used for integral anti-windup and dispatch)
            %   - vBus [V], is the voltage of each bus in the medium voltage grid
            %
            % Output:
            %   - u [-], is the set-points for the controllable assets
            
            % Control system (proportional control with feedforward)
            % Check if sampling time slower than asset dynamics
            if obj.Ts<obj.tauMin
                xc = obj.Ctrld.a*obj.xcOld + obj.Ctrld.b*[-(zhat+ref);sum(yProd)];
                uP = obj.Ctrld.c*obj.xcOld + obj.Ctrld.d*[-(zhat+ref);sum(yProd)];
                obj.xcOld = xc;
            else
                uP = obj.kP*(1/obj.numProd)*ones(obj.numProd,1)*(-(zhat+ref));
            end
 
            % Integral control
            if obj.wInt == true
                % Calculate error
                err = ref-yMeas;
                % Calculate control signal
                u = uP + ones(obj.numProd,1)*obj.uIOld;
                % Antiwindup
                for i=1:obj.numProd
                    if u(i) > pMax(i)
                        u(i) = pMax(i);
                        obj.satFlagUp(i) = 1;
                    else
                        obj.satFlagUp(i) = 0;
                    end

                end
                for i=1:obj.numProd
                    if u(i) < pMin(i)
                        u(i) = pMin(i);
                        obj.satFlagDown(i) = 1;
                    else
                        obj.satFlagDown(i) = 0;
                    end
                end
                % Check saturation flags
                if sum(obj.satFlagUp) >= obj.numProd
                    obj.uIOld = obj.uIOld;
                elseif sum(obj.satFlagDown) >= obj.numProd
                    obj.uIOld = obj.uIOld;
                else
                    % Integral control
                    obj.uIOld = obj.uIOld + obj.kI*err; 
                end
                % Calculate new dispatch
                if obj.wDispatch
                    pRef = sum(u);%+obj.uIOld;
                    uDis = dispatch(obj,pRef,pMax,pMin,yhat,vBus);
                    u = uDis;
                end
            % Without integral control    
            else
                % Calculate dispatch
                if obj.wDispatch
                    pRef = sum(uP);
                    uDis = dispatch(obj,pRef,pMax,pMin,yhat,vBus);
                    u = uDis;
                % Without Dispatch    
                else
                    u = uP;
                end
            end
            obj.uOld = u;
        end
        
        % Dispatch
        function uDis = dispatch(obj,pRef,pMax,pMin,pLoad,vBus)
            % The dispatch algorithm uses assumptions of knowledge of power
            % consumption at load busses, measurements of voltages. The
            % dispatch optimization problem is then convex and cvx is used
            % to solve the problem.
            %
            % Input:
            %   - pRef [W], is the power reference from the feedforward and
            %     feed back control.
            %   - pMax [W], is a vector of maximum allowable power of
            %     flexible assets.
            %   - pMin [W], is a vector of minimum allowable power of
            %     flexible assets.
            %   - pLoad [W], is a vector of estimated power of inflexible
            %     loads.
            %   - vBus [V], is a vector of measured voltage on all busses.
            %
            % Output:
            %   - uDis [W], is a vector of dispatched active power
            %     references for the flexible assets.
                        
            % Setup optimization problem 
            sf = max(abs(pLoad));       % Numerical problems so needs scaling
            v = vBus/sf;                % Scale voltages at each bus
            V1 = diag(1./v);            % Form reciprocal voltage matrix
            B = real(V1*obj.Z*V1);      % Calculate real part of loss matrix
            pRef = pRef/sf;             % Scale the reference signal
            pMax = pMax/sf;             % Scale the maximum values
            pMin = pMin/sf;             % Scale the minimum values
            pLoad = pLoad/sf;           % Scale the load values
            P = zeros(obj.numBus,1);    % Allocate the power vector
            % Identify index of each bus (0 = free, 1 = load, 2 = flexible)
            idxNotProd = [obj.idFree obj.idLoad];
            idxLoad = obj.idLoad;
            idxProd = obj.idProd;
            % Input estimated load values into power vector
            P(idxLoad) = pLoad;

            % Weights
            w1 = 0.01;     % Active power loss
            w2 = 2000;    % Reference following
            w3 = 0.5;    % Rate change
            
            % Solve optimization problem using CVX (see http://cvxr.com for more information)
            cvx_begin quiet
                % Set variables
                variables Pp(obj.numBus)
                % Cost function with weights
                minimize(w1*(Pp'*B*Pp+2*P'*B*Pp + P'*B*P) + w2*norm(sum(Pp)-pRef,2) + w3*norm(obj.PpOld-Pp(idxProd),1))
                % Constraints
                subject to
                    % All inflexible busses must be zero
                    Pp(idxNotProd) == 0;
                    % Flexible busses must be within constraints
                    pMin<=Pp(idxProd)<=pMax;
            cvx_end
            % Save output, used for rate minimization
            obj.PpOld = Pp(idxProd);
            % Set output
            uDis = Pp(idxProd)*sf;
        end
        
        

        %% Set and get functions
        % Set the proportional gain
        function setProGain(obj,kP)
            obj.kP = kP;
        end
        % Set the integral gain
        function setIntGain(obj,kI)
            obj.kI = kI;
        end
        % Tuning of the Kalman filter in the observer
        function tuneDistEKF(obj,Q,R)
            obj.Qd = Q;     % State covariance
            obj.Rd = R;     % Output covariance
        end
        
        % Turn on/off integral control
        function withIntegral(obj,in)
            % Input:
            %   - in [-], true or false
            obj.wInt = in;
        end
        
        % Turn on/off optimal dispatch
        function withDispatch(obj,in)
            % Input:
            %   - in [-], true of false
            obj.wDispatch = in;
        end
           
    end
    
end

% Internal functions
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
At = zeros(2,2);
for i = 1:numOsci
    A = blkdiag(A,At);
end
% Finalyse system matrix by inserting for bias term
A = blkdiag(A,0);
A = expm(A*Ts); % Discritise
% Sample system
x = A*xold;

% Form output matrix
C = zeros(1,numOsci*4+1);
for i=1:numOsci*2
    C(i) = x(numOsci*2+i);
end
C(end) = 1;

y = C*xold;
end

