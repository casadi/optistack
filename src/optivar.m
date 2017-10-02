classdef optivar < OptimizationObject
    % Create an OPTISTACK optimization variable
    %
    %  x = optivar()                   Scalar
    %  x = optivar(n)                  Column vector of length n
    %  x = optivar(n,m)                Matrix of shape n-by-m
    %  x = optivar(n,m,name)           Supply a name for printing
    %  x = optivar(MX mx)              Wrap an existing MX object
    %
    %  optivar Methods:
    %    optivar  - constructor
    %    setInit  - Set initial value for use in NLP solver
    %    setLb    - Give a lower bound to the optimization variable
    %    setUb    - Give an upper bound to the optimization variable
    
    properties(Constant, Hidden=true)
       shorthand = 'x'; 
    end
    
    properties(SetAccess=private, Hidden=true)
       lb
       ub
       init
    end

    methods
        function [] = setInit(self,v)
            % Set initial value for use in NLP solver
            %
            % Supplied value v may be either a scalar or a matrix with
            % matching shape
            self.init(:,:) = v;
        end
        function [] = setLb(self,v)
            % Give a lower bound to the optimization variable
            %
            % Supplied value v may be either a scalar or a matrix with
            % matching shape
            self.lb(:,:) = v;
        end
        function [] = setUb(self,v)
            % Give an upper bound to the optimization variable
            %
            % Supplied value v may be either a scalar or a matrix with
            % matching shape
            self.ub(:,:) = v;
        end
        function self = optivar(varargin)
            % Create an OPTISTACK optimization variable
            %  x = optivar()                   Scalar
            %  x = optivar(n)                  Column vector of length n
            %  x = optivar(n,m)                Matrix of shape n-by-m
            %  x = optivar(n,m,name)           Supply a name for printing
            %  x = optivar(obj)                Wrap an existing MX object
           if nargin == 0
               shape = 1;
           elseif nargin == 1
               shape = varargin{1};
           elseif nargin > 1
               shape = [varargin{1};varargin{2}];
           end
           
           if nargin == 3
               name = varargin(3);
           else
               if isnumeric(shape)
                   name = {'x'};    % automatic name
               else
                   name = {};       % wrapping an MX object
               end
           end
           
           import casadi.*
           self@OptimizationObject(shape,name{:});
           self.lb = -inf*DM.ones(self.sparsity());
           self.ub = inf*DM.ones(self.sparsity());
           self.init = DM.zeros(self.sparsity());
        end
    end
    
end

