
classdef OptimizationObject < casadi.MX
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        mapping = containers.Map('KeyType','uint64','ValueType','any');
    end
    properties
        value = nan;
    end
    
    methods
        function n = numArgumentsFromSubscript(self,s,callingContext)
           n=1;
        end
        function self = OptimizationObject(shape,name)
           if isscalar(shape)
               shape = [shape 1];
           end
           self@casadi.MX(casadi.MX.sym(name,casadi.Sparsity.dense(shape(1),shape(2))));
           mymapping = OptimizationObject.mapping;
           mymapping(self.hash()) = self;
        end
        function val = optival(self)
            val = self.value;
        end
        function [] = setValue(self,value)
            self.value = value;
        end
        function r = reshape(self,a,b)
           a
           b
           if size(self,1) == a && size(self,2)==b
              r = self;
           else
              r = casadi.reshape(self,a,b);
           end
        end
    end
    methods(Static)
        function [ syms ] = get_primitives( el, varargin )
            % Out of a list of expression, retrieves all primitive expressions
            % The result is sorted into a dictionary with the key originating
            % from the 'shorthand' class attribute of OptimzationObject subclasses
            if isempty(varargin)
               dep = true;
            else
               dep = varargin{1}; 
            end

            orig_state = warning;
            warning('off','all')

            try
                vars = symvar(veccat(el{:}));
            catch

                vars = {}; 
            end
            warning(orig_state);
            
            syms = struct();

            for i=1:length(vars)
                v = vars{i};
                mymapping = OptimizationObject.mapping;
                vobj = mymapping(v.hash());
                symkey = vobj.shorthand;
                if isfield(syms,symkey)
                    syms.(symkey) = {syms.(symkey){:}, vobj};
                else
                    syms.(symkey) = { vobj };
                end
            end
        end
    end  
end

