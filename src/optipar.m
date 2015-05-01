classdef optipar < OptimizationObject
    % Create an optimization parameter,
    % something that is not optimised for.
    %
    %  x = optipar()                   Scalar
    %  x = optipar(n)                  Column vector of length n
    %  x = optipar(n,m)                Matrix of shape nxn
    %  x = optipar(n,m,name)           A name for printing the variable
    
    properties(Constant)
       shorthand = 'p';
    end

    methods
        function self = optipar(varargin)
           if isempty(varargin)
               shape = 1;
           elseif length(varargin)==1
               shape = varargin{1};
           elseif length(varargin)==2
               shape = [varargin{1};varargin{2}];
           end
           
           if length(varargin)==3
               name = varargin{3};
           else
               name = 'p';
           end
           
           self@OptimizationObject(shape,name);
        end
    end
    
end

