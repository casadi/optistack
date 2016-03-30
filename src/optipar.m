classdef optipar < OptimizationObject
    % Create an OPTISTACK optimization parameter,
    % something that is not optimized for.
    %
    %  x = optipar()                   Scalar
    %  x = optipar(n)                  Column vector of length n
    %  x = optipar(n,m)                Matrix of shape n-by-m
    %  x = optipar(n,m,name)           A name for printing the variable
    %
    %  optipar Methods:
    %    setValue  - Set the parameter value for use in NLP solver
    properties(Constant, Hidden=true)
       shorthand = 'p';
    end

    methods
        function [] = setValue(self,v)
            % Set the parameter value for use in NLP solver
            %
            % Supplied value v may be either a scalar or a matrix with
            % matching shape
            setValue@OptimizationObject(self,v);
        end
        function self = optipar(varargin)
            % Create an OPTISTACK optimization parameter,
            % something that is not optimized for.
            %
            %  x = optipar()                   Scalar
            %  x = optipar(n)                  Column vector of length n
            %  x = optipar(n,m)                Matrix of shape n-by-m
            %  x = optipar(n,m,name)           A name for printing the variable
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
               name = 'p';
           end
           
           self@OptimizationObject(shape,name);
        end
    end
    
end

