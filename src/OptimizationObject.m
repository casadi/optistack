classdef OptimizationObject < casadi.MX
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function self = OptimizationObject(shape,name)
           if isscalar(shape)
               shape = [shape 1];
           end
           self@casadi.MX(casadi.MX.sym(name,casadi.Sparsity.dense(shape(1),shape(2))));
        end
    end
    
end

