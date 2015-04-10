classdef optivar < casadi.MX
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function self = optivar(foo)
           self@casadi.MX(casadi.MX.sym('foo'));
        end
    end
    
end

