classdef lvgcM2 < handle
    %LVGCM2 (low voltage grid controller - Mark 2) is a low voltage
    % controller class for the simulation framework DiSC.
    %
    %   This class will implement a specific controller class for a low
    %   voltage grid. Its assignments are to make a LV grid follow a power
    %   reference and apply voltage control.
    %   The controller is based on a consensus algorithm, which
    %   uses the systems communication graph. Thus, the key variables are
    %   the communication Laplacian matrix, L, and the weights alpha.
    %
    % Revision:
    % 06-05-2015, R. Pedersen, Aalborg University. Notes:

    % Parameters for the controller object
    properties
        % General parameters
        vBase = 400;            % [V]. Voltage base.
        numAssets = 0;          % [-]. Number of flexible units.
        Ts = 60;                % [s]. Sampling time.
        
        % Reference following control parameters
        alpha = [1 0.5 0.25];   % [-]. Laplacian weights.
        L                       % [-]. Communication Laplacian graph.
        K = 1;                  % [-]. Proportional gain
        Ki = 0;                 % [-]. Integral gain
        intErr = 0;             % [-]. Integral error
        weightPos = 1;          % [-]. Flag indicating which value in alpha should be used {1,...,n}. Default is {1,2,3}.
        oldWeightPos = 1;       % [-]. Former weight for the Laplacian, used for guard on control switching.
        weightSwitchTimer = 0;  % [-]. Timer value for determining the last time a switch has occured. 
        pAssetMax               % [W]. Vector of maximum values of assets.
        pAssetMin               % [W]. Vector of minimum values of assets.
        satFlagUp               % [-]. Saturation flags used for anti-windup
        satFlagDown             % [-]. Saturation flags used for anti-windup
        hMaxStart = 1.06;       % [PU]. Hysterisis for voltage control start.
        hMaxStop = 1.05;        % [PU]. Hysterisis for voltage control stop.
        hMinStart = 0.94;       % [PU]. Hysterisis for voltage control start.
        hMinStop = 0.95;        % [PU]. Hysterisis for voltage control stop.
        vFlag                   % [-]. Flag used in the voltage control.
        vDroopGains             % [-]. Vector of droop gains used for voltage control.
        aOld = 1;               % [-]. Initial laplacian weight.
        pAssetOld               % [-]. Former asset power values, used for anti-windup in case of com. loss.
    end
    
    methods
        % Constructor
        function obj = lvgcM2(param)
            % param format:
            % - param.numAssets     [-]. Number of assets the LVGC can control
            % - param.L             [-]. Communication Laplacian
            % - param.Ts            [s]. Sampling time
            
            % Set general parameters
            obj.vBase = param.vBase;
            obj.numAssets = param.numAssets;
            obj.L = param.L;
            obj.Ts = param.Ts;
            obj.pAssetMax = param.pAssetMax;
            obj.pAssetMin = param.pAssetMin;
            obj.vDroopGains = param.vDroopGains;
            
            obj.vFlag = zeros(param.numAssets,1);
            obj.satFlagUp = zeros(obj.numAssets,1);
            obj.satFlagDown = zeros(obj.numAssets,1);
            obj.pAssetOld = zeros(obj.numAssets,1);
        end
        
        % Sample the controller
        function [uRefP,uRefQ,pFlexUpAgg,pFlexDownAgg,qFlexUpAgg,qFlexDownAgg] = sample(obj,pRef,pMeas,pAsset,pMin,pMax,qMin,qMax,vBus)
            
            % Power reference control
            [uRefP,pFlexUpAgg,pFlexDownAgg] = pRefFollow(obj,pRef,pMeas,pAsset,pMin,pMax);
            
            % Voltage Control
            [uRefQ,qFlexUpAgg,qFlexDownAgg] = voltControl(obj,vBus,qMin,qMax);
        end

        % Power reference following controller
        function [uRefP,pFlexUpAgg,pFlexDownAgg] = pRefFollow(obj,pRef,pMeas,pAsset,pMin,pMax)
            % Error
            err = pRef-pMeas;
            % Set weight on Laplacian matrix
            a = weightLaplacian(obj);
            if a~=obj.aOld
                temp1 = (obj.Ki*obj.aOld)*obj.intErr;
                obj.intErr = temp1/(obj.Ki*a);
            end
            obj.aOld = a;
            
            % Control output
            pOut = obj.K*err + obj.Ki*obj.intErr;
            % Divide between assets
            if pOut<=0
                pOut = pOut/abs(sum(obj.pAssetMin));
                x = [pOut;pAsset./abs(obj.pAssetMin)];
            else
                pOut = pOut/abs(sum(obj.pAssetMax));
                x = [pOut;pAsset./abs(obj.pAssetMax)];
            end
            % Apply "consensus"        
 
            x1 = x-a*(obj.L*x);
            if pOut <= 0
                uRefP = x1(2:end).*abs(obj.pAssetMin);
            else
                uRefP = x1(2:end).*abs(obj.pAssetMax);
            end
            
            % Anti-windup
            for i=1:obj.numAssets
                if pAsset(i)==obj.pAssetOld(i)
                        obj.satFlagUp(i) = 1;
                        obj.satFlagDown(i) = 1;
                        if uRefP(i) > pMax(i)
                            uRefP(i) = pMax(i);
                        elseif uRefP(i) < pMin(i)
                            uRefP(i) = pMin(i);
                        end
                else
                    if uRefP(i) > pMax(i)
                        uRefP(i) = pMax(i);
                        obj.satFlagUp(i) = 1;
                    elseif uRefP(i) < pMin(i)
                        uRefP(i) = pMin(i);
                        obj.satFlagDown(i) = 1;
                    else
                        obj.satFlagUp(i) = 0;
                        obj.satFlagDown(i) = 0;
                    end
                end
            end
            % Check saturation flags
            if sum(obj.satFlagUp) >= obj.numAssets || sum(obj.satFlagDown) >= obj.numAssets
                if err >= 0 && obj.intErr <= 0
                    obj.intErr = obj.intErr + err;
                elseif err < 0 && obj.intErr > 0
                    obj.intErr = obj.intErr + err;
                else
                    obj.intErr = obj.intErr;
                end
            else
                % Integral control
                obj.intErr = obj.intErr + err;
            end
            
            % Aggregate Flexibility
            pFlexUpAgg = sum(pMax)-sum(uRefP);
            pFlexDownAgg = sum(pMin)-sum(uRefP);
        end
        
        % Voltage control
        function [uRefQ,qFlexUpAgg,qFlexDownAgg] = voltControl(obj,vMeas,qFlexUp,qFlexDown)
            uRefQ = zeros(obj.numAssets,1);
            % Hysteresis for starting voltage control (0~ok, 1~high voltage, -1~low voltage)
            vFlagTemp = ((obj.vFlag == 1) & ((abs(vMeas)'/obj.vBase)>obj.hMaxStop)) | ((abs(vMeas)'/obj.vBase)>obj.hMaxStart);
            obj.vFlag = vFlagTemp - (((obj.vFlag == -1) & ((abs(vMeas)'/obj.vBase)<obj.hMinStop)) | ((abs(vMeas)'/obj.vBase)<obj.hMinStart));
            
            for i=1:obj.numAssets
                if obj.vFlag(i)~=0
                    uRefQ(i) = obj.vDroopGains(i)*(1-abs(vMeas(i))/obj.vBase);
                end
            end
            
            % Aggregate Flexibility
            qFlexUpAgg = sum(qFlexUp);
            qFlexDownAgg = sum(qFlexDown);
        end
        
        % Determine which weight for the Laplacian should be used
        function a = weightLaplacian(obj)
            switchTime = obj.Ts*5;  % Delay 5 samples
            obj.weightSwitchTimer = obj.weightSwitchTimer + obj.Ts;
            % Check if the weight has changed
            if obj.weightPos > obj.oldWeightPos
                a = obj.alpha(obj.weightPos);
                obj.oldWeightPos = obj.weightPos;
                obj.weightSwitchTimer = 0;
            elseif obj.weightPos < obj.oldWeightPos && obj.weightSwitchTimer > switchTime
                a = obj.alpha(obj.oldWeightPos-1);
                obj.oldWeightPos = obj.oldWeightPos-1;
                obj.weightSwitchTimer = 0;
            else
                a = obj.alpha(obj.oldWeightPos);
            end
        end
        
        %% Set and get functions
        % Set the proportional and integral gains
        function setGains(obj,K,Ki)
            obj.K = K;
            obj.Ki = Ki;
        end
        
        % Set communication Laplacian 
        function setLaplacian(obj,L)
            % NOTE: that for now the Laplacian must have a certain
            % structure. The LVGC is communicating with all assets, and
            % there is no asset - asset communication.
            obj.L = L;
        end
        
        % Set alpha gain
        function setLaplacianWeight(obj,alpha)
            % Input:
            % - alpha [-], is a set of Laplacian weights, [a1 a2 ... an]. Three as default
            obj.alpha = alpha;
        end
        
        % Set which value in the alpha vector should be used as weight
        function setLapWeightPos(obj,in)
            in = fix(in); % Ensures that the position is read as an integer.
            if in <=0 || in > length(obj.alpha)
                error('in value needs to be an integer between 1 and n, where n is the number of values in the set of Laplacian weights');
            end
            obj.weightPos = in;
        end
        
        % Set hysterisis on voltage control start/stop
        function setHystParam(obj,hMaxStart,hMaxStop,hMinStart,hMinStop)
            obj.hMaxStart = hMaxStart;
            obj.hMaxStop = hMaxStop;
            obj.hMinStart = hMinStart;
            obj.hMinStop = hMinStop;
        end
           
    end
    
end
