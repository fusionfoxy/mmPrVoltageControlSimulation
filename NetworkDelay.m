classdef NetworkDelay
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        distribution;
        Ts;
    end
    
    methods
        function obj = NetworkDelay(param)
            obj.distribution = param.mu;
            obj.Ts = param.Ts;
        end
        
        function d = draw(obj)
            d = exprnd(obj.distribution);
            d = ceil(d/obj.Ts);
        end
    end
    
end

