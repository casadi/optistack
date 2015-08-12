classdef optivar < OptimizationObject
    % Create an optimization variable
    %
    %  x = optivar()                   Scalar
    %  x = optivar(n)                  Column vector of length n
    %  x = optivar(n,m)                Matrix of shape nxn
    %  x = optivar(n,m,name)           A name for printing the variable
    
    properties(Constant)
       shorthand = 'x'; 
    end
    
    properties
       lb = -inf;
       ub = inf;
       init = 0;
    end

    methods
        function [] = setInit(self,v)
            self.init = v;
        end
        function [] = setLb(self,v)
            self.lb = v;
        end
        function [] = setUb(self,v)
            self.ub = v;
        end
        function self = optivar(varargin)
           if isempty(varargin)
               shape = 1;
           elseif length(varargin)==1
               shape = varargin{1};
           elseif length(varargin)>1
               shape = [varargin{1};varargin{2}];
           end
           
           if length(varargin)==3
               name = varargin{3};
           else
               name = 'x';
           end
           
           self@OptimizationObject(shape,name);
        end
    end
    
end

