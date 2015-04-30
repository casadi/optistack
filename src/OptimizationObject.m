
classdef OptimizationObject < casadi.MX
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        mapping = containers.Map('KeyType','uint64','ValueType','any');
    end
    properties
        num_value = nan;
    end
    
    methods
        function self = OptimizationObject(shape,name)
           if isscalar(shape)
               shape = [shape 1];
           end
           self@casadi.MX(casadi.MX.sym(name,casadi.Sparsity.dense(shape(1),shape(2))));
           mymapping = OptimizationObject.mapping;
           mymapping(self.hash()) = self;
        end
        function val = optival(self)
            val = self.num_value;
        end
        function [] = setValue(self,value)
            self.num_value = value;
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

            vars = getSymbols(el{:});
            
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

