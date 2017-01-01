classdef lvgcM1 < handle
    %lvgcM1 is a low voltage controller class for the simulation
    %framework DiSC.
    %   This class will implement a specific controller class for a low
    %   voltage grid. Its assignments are to ensure voltages at each bus
    %   and make the entire LV grid follow an energy reference. Following
    %   an energy reference can be used for e.g., engaging in the energy
    %   markets or offering frequency control services. In practise the
    %   controller is designed to follow a power reference.
    %
    %   The controller ss based on a "fairness" consensus algorithm which
    %   can be applied when the communication graph structure is known. The
    %   "fairness" stems from the controllable assests in the grid equally
    %   contribute compared to their size.
    
    properties
        vBase = 400;        % [V]. Voltage base.
        kP = 0.3;           % [-]. Proportional gain.
        kI = 0.2;           % [-]. Integral gain.
        errI = 0;           % [-]. Integral error.
        numAssets = 0;      % [-]. Number of flexible assets
        uRefOld             % [-]. Plaveholdet for output
        hMaxStart = 1.06;   % [PU]. Hysterisis for voltage control start
        hMaxStop = 1.05;    % [PU]. Hysterisis for voltage control stop
        hMinStart = 0.94;   % [PU]. Hysterisis for voltage control start
        hMinStop = 0.95;    % [PU]. Hysterisis for voltage control stop
        vFlag               % [-]. Flag used in the voltage contro
        vDroopGains
        pFlexUpOld
        pFlexDownOld
        qFlexUpOld
        qFlexDownOld
        
    end
    
    methods
        % Constructor
        function obj = lvgcM1(param)
            % param format:
            % - vBase
            % - numAssets
            obj.vBase = param.vBase;
            obj.numAssets = param.numAssets;
            obj.vDroopGains = param.vDroopGains;
            obj.uRefOld = zeros(param.numAssets,1);
            obj.vFlag = zeros(param.numAssets,1);
            obj.pFlexUpOld = zeros(param.numAssets,1);
            obj.pFlexDownOld = zeros(param.numAssets,1);
            obj.qFlexUpOld = zeros(param.numAssets,1);
            obj.qFlexDownOld = zeros(param.numAssets,1);
        end
        
        % Sample the controller
        function [uRefP,uRefQ,pFlexUpAgg,pFlexDownAgg,qFlexUpAgg,qFlexDownAgg] = sample(obj,pRef,pMeas,pFlexUp,pFlexDown,qFlexUp,qFlexDown,satFlagUp,satFlagDown,vMeas)
            
            uRefP = refFollow(obj,pRef,pMeas,satFlagUp,satFlagDown);
            uRefQ = voltCtrl(obj,vMeas);
            [pFlexUpAgg,pFlexDownAgg,qFlexUpAgg,qFlexDownAgg] = aggFlex(obj,pFlexUp,pFlexDown,qFlexUp,qFlexDown);
        end
        
        % Reference following
        function uRefP = refFollow(obj,pRef,pMeas,satFlagUp,satFlagDown)
            uRefP = zeros(obj.numAssets,1);
            % Calculate error
            err = pRef-pMeas;
            
            % Check if input measurement is NaN (this indicates lost package)
            satUpTemp= satFlagUp;
            satDownTemp  = satFlagDown;
            nanIdx = isnan(satFlagUp);
            satUpTemp(nanIdx) = 1;          % If package is lost assume saturated
            nanIdx = isnan(satFlagDown);
            satDownTemp(nanIdx) = 1;        % If package is lost assume saturated
            
            % Check for asset availability
            numAvabAssets = obj.numAssets - sum(satUpTemp) + sum(satDownTemp);
            % Anti-windup
            if numAvabAssets ~= 0
                obj.errI = obj.errI + obj.kI*err;
            else
                obj.errI = obj.errI-sign(obj.errI)*abs(err);
            end
            % Set control output
            for i=1:obj.numAssets
                if isnan(satFlagUp(i)) || isnan(satFlagDown(i))
                    uRefP(i) = obj.uRefOld(i);
                else
                    uRefP(i) = (obj.kP*err + obj.errI)/obj.numAssets;
                end
                
            end
%             for i=1:obj.numAssets
%                 if satFlagUp(i)==0 && satFlagDown(i)==0
%                     uRefP(i) = (obj.kP*err + obj.errI)/numAvabAssets;
%                 else
%                     uRefP(i) = obj.uRefOld(i);
%                 end
%             end
            obj.uRefOld = uRefP;
        end
        
        % Voltage control
        function uRefQ = voltCtrl(obj,vMeas)
            uRefQ = zeros(obj.numAssets,1);
            % Hysteresis for starting voltage control (0~ok, 1~high voltage, -1~low voltage)
            vFlagTemp = ((obj.vFlag == 1) & ((abs(vMeas)'/obj.vBase)>obj.hMaxStop)) | ((abs(vMeas)'/obj.vBase)>obj.hMaxStart);
            obj.vFlag = vFlagTemp - (((obj.vFlag == -1) & ((abs(vMeas)'/obj.vBase)<obj.hMinStop)) | ((abs(vMeas)'/obj.vBase)<obj.hMinStart));
            
            for i=1:obj.numAssets
                if obj.vFlag(i)~=0
                    uRefQ(i) = obj.vDroopGains(i)*(1-abs(vMeas(i))/obj.vBase);
                end
            end

        end
        
        % Aggregate flexibility
        function [pFlexUpAgg,pFlexDownAgg,qFlexUpAgg,qFlexDownAgg] = aggFlex(obj,pFlexUp,pFlexDown,qFlexUp,qFlexDown)
            
            for i=1:obj.numAssets
                if isnan(pFlexUp(i))
                    pFlexUp(i) = obj.pFlexUpOld(i);
                end
                if isnan(pFlexDown(i))
                    pFlexDown(i) = obj.pFlexDownOld(i);
                end
                if isnan(qFlexUp(i))
                    qFlexUp(i) = obj.qFlexUpOld(i);
                end
                if isnan(qFlexDown(i))
                    qFlexDown(i) = obj.qFlexDownOld(i);
                end       
            end
            pFlexUpAgg = sum(pFlexUp);
            pFlexDownAgg = sum(pFlexDown);
            qFlexUpAgg = sum(qFlexUp);
            qFlexDownAgg = sum(qFlexDown);
            
            obj.pFlexUpOld = pFlexUp;
            obj.pFlexDownOld = pFlexDown;
            obj.qFlexUpOld = qFlexUp;
            obj.qFlexDownOld = qFlexDown;
            
        end
        
        % Set and Get functions
        function setGains(obj,kP,kI)
            obj.kP = kP;
            obj.kI = kI;
        end
        
    end
end