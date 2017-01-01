classdef mvgcM2 < handle
    %MVGCM2 (medium voltage grid controller - Mark 2) is a medium voltage
    % controller class for the simulation framework DiSC.
    %
    %   This class will implement a specific controller class for a medium
    %   voltage grid. Its assignments are to ensure voltages at each bus
    %   and make the entire MV grid follow an energy reference. Following
    %   an energy reference can be used for e.g., engaging in the energy
    %   markets or offering frequency control services. In practise the
    %   controller is designed to follow a power reference. To implement
    %   the controller Model Predictive Control (MPC) is used.
    %
    %   The controller also incorporates an observer for estimating the
    %   power consumption of uncontrollable loads and production from
    %   controllable production units. Controllable production units could
    %   be a wind or solar power plant. To estimate uncontrollable 
    %   aggregated consumption, oscillator models are used.
    %
    %   Flexible units are modelled by first order systems:
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
        numAsset = 0;       % [-]. Number of flexible assets.
        numLoad = 0;        % [-]. Number of aggregated loads.
        numOsci = 2;        % [-]. Number of oscillators for each load model.
        Ts = 60;            % [s]. Sampling time.
        busType             % [-]. Type of each bus in the system (0 = Free, 1 = Load, 2 = Asset).
        
        % Model parameters
        % Flexible assets
        sysc_d              % [-]. Flexible asset system (discrete time).
        
        % Observer parameters
        Qd                  % [-]. State covariance of disturbance, can be used for tuning.
        Rd                  % [-]. Output covariance matrix, ca be used for tuning.
        PdOld               % [-]. Old state covariance, used for the Kalman filter.
        xhatOld             % [-]. Holds the old estimates of states.
        numDstates          % [-]. Number of disturbance states.
        numCstates          % [-]. Number of controllable states.
        
        % Control parameters
        H = 60*60*6         % [-]. MPC horizon
        tauMin = 10;        % [sec]. Minimum time constant for applying dynamic control.
        kI = 0.08;           % [-]. Integral gain
        xInt = 0;
        
        % Dispatch parameters
        Z                   % [-]. Is the systems bus impedance matrix, used for loss minimization.
        numBus = 0          % [-]. Number of busses in the system.
        idFree              % [-]. Index of free busses (no load or asset).
        idLoad              % [-]. Index of load busses.
        idAsset             % [-]. Index of flexible asset busses.
    end
    
    methods
        % Constructor
        function obj = mvgcM2(param)
            % param format:
            % - param.aProd         [-]. Is a vector containing the 'a' parameter of flexible assets
            % - param.bProd         [-]. Is a vector containing the 'b' parameter of flexible assets
            % - param.cProd         [-]. Is a vector containing the 'c' parameter of flexible assets
            % - param.numLoad       [-]. Number of uncontrollable aggregated loads
            % - param.Ts            [s]. Sampling time
            
            % Set general parameters
            obj.numAsset = length(param.aAsset);
            obj.numLoad = param.numLoad;
            obj.Z = param.Z;
            obj.Ts = param.Ts;
            obj.busType = param.busType;
            obj.numDstates = obj.numLoad*(obj.numOsci*4+1);
            obj.numCstates = obj.numAsset;
            
            % Set parameters for observer
            % Disturbance estimator
            obj.Qd = eye(obj.numDstates);                   % State covariance (Tuning parameter)
            obj.Rd = eye(obj.numLoad);                      % Output covariance (Tuning parameter)
            obj.PdOld = zeros(obj.numDstates);              % Kalman filter covariance
            obj.xhatOld = ones(obj.numCstates+obj.numDstates,1);
            
            % Set parameters controller
            obj.tauMin = min(abs(1./param.aAsset))/15;
            Ac = diag(param.aAsset);
            Bc = diag(param.bAsset);
            Cc = diag(param.cAsset);
            Dc = zeros(size(Cc,1),size(Bc,2));
            sysc = ss(Ac,Bc,Cc,Dc);
            obj.sysc_d = c2d(sysc,obj.Ts);
            obj.H = 5*60*60/obj.Ts;
            
            % Set parameters for dispatch
            obj.numBus = length(obj.busType);
            obj.idFree = find(obj.busType==0);
            obj.idLoad = find(obj.busType==1);
            obj.idAsset = find(obj.busType==2);
        end
        
        % Sample the controller
        function [yhat,xhat,u]= sample(obj,Pref,yAsset,yDist,yMeas,Pmax,Pmin,vBus)
                      
            % Estimation
            [xhat,yhat,zhat,Ad,Cd] = observer(obj,yDist);
            
             % Control
             u = controller(obj,Ad,Cd,Pref,yAsset,yMeas,xhat,Pmax,Pmin,vBus);
        end
        
        % Observer
        function [xhat,yhat,zhat,Ad,Cd] = observer(obj,y)
            % The observer implements an extended Kalman filter for
            % estimating the states of the controllable system and the
            % states of the uncontrollable disturbance. Also, it estimates
            % some of the parameters of the disturbance model. This is way
            % it is implemented as an extended Kalman filter.
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
        function u = controller(obj,Ad,Cd,Pref,yAsset,yMeas,xhat,Pmax,Pmin,vBus)
            % The controller is implemented as an internal model controller
            % which uses the disturbance model in the control apporach.
            % Further, the system has been augmented with an additional
            % integral state. The integral controller uses the measurements
            % from the controllable part directly, to not have any bias
            % which may occure if only using the state estimates. 
            %
            % Input:
            %   - Ad [-], is the disturbance model A matrix
            %   - Cd [-], is the disturbance model C matrix
            %   - Pref [W], is the active power reference signal for the controller
            %   - yAsset [W], is a vector of active power measurements at the flexible assets
            %   - yMeas [W], is the measurement of the entire systems power (used for integral control)
            %   - xhat [W], is a vector containing the estimated disturbance state
            %   - Pmax [W], is a vector containing information on the maximum power available
            %   - Pmin [W], is a vector containing information on minimum power available
            %   - vBus [V], is the voltage of each bus in the medium voltage grid
            %
            % Output:
            %   - u [-], is the set-points for the controllable parts
            
            % Setup MPC optimization problem
            Ac = obj.sysc_d.A;
            Bc = obj.sysc_d.B;
            Cc = obj.sysc_d.C;
            % Integral action
            obj.xInt = (Pref(1)-yMeas);
            % Solve MPC optimization problem
            cvx_begin quiet
                % Define decision variables 
                variables xd(obj.numDstates,obj.H) xc(obj.numAsset,obj.H) uP(obj.numAsset,obj.H-1)
                % Objective function
%                 minimize(norm((Pref)-(sum(Cc*xc(:,1:end))-sum(Cd*xd(:,1:end))),2))
                minimize(norm((Pref+obj.xInt)-(sum(Cc*xc(:,1:end))),2))
                % Constraints
                subject to
                    % Disturbance model
                    xd == [xhat Ad*xd(:,1:end-1)];
                    
                    % Control model
                    if obj.Ts<obj.tauMin
                        xc == [yAsset Ac*xc(:,1:end-1)+Bc*uP];
                    else
                        xc == [yAsset uP]; 
                    end
                    
                    % Control output constraints
                    repmat(Pmin,1,obj.H-1) <= uP <= repmat(Pmax,1,obj.H-1);    
            cvx_end
            % Dispatch output
            u = uP(:,1);
            
        end
        
        % Dispatch
        function uDis = dispatch(obj,Pref,Pmax,Pmin,Pload,vBus)
            sf = max(abs(Pload));
            V = diag(abs(1./vBus));
            R = real(obj.Z);
            B = V*R*V;
            Pref = Pref/sf;
            Pmax = Pmax/sf;
            Pmin = Pmin/sf;
            Pload = Pload/sf;
            P = zeros(obj.numBus,1);
            idxNotProd = [obj.idFree obj.idLoad];
            idxLoad = [obj.idLoad obj.idBoth];
            idxProd = [obj.idProd obj.idBoth];
            P(idxLoad) = Pload;
            % Form load currents
            cvx_begin quiet
                variables Pp(obj.numBus)
                %minimize(Pp'*B*Pp+2*P'*B*Pp + P'*B*P + norm(obj.PpOld-Pp) + 1.2*abs(Pp(4)) + 1.2*abs(Pp(3)))
                minimize(Pp'*B*Pp+2*P'*B*Pp + P'*B*P + 6*norm(sum(Pp)-Pref,2) + norm(obj.PpOld-Pp(idxProd),1) + [0.1 0.1 0.1]*abs(Pp(idxProd)))
                subject to
                    Pp(idxNotProd) == 0;
                    %sum(Pp) == Pref;
                    Pmin<=Pp(idxProd)<=Pmax;
            cvx_end
            Pp(isnan(Pp)) = 0;
            obj.PpOld = Pp(idxProd);
            uDis = Pp(idxProd)*sf;
        end
        
        

        %% Set and get functions
        % Tuning of the Kalman filter of the observer
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

