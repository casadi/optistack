
classdef OptimizationObject < casadi.MX
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant, Hidden=true)
        mapping = containers.Map('KeyType','uint64','ValueType','any');
    end
    properties(Hidden=true)
        value = nan;
    end
    
    methods(Hidden=true)
        function n = numArgumentsFromSubscript(self,s,callingContext)
           n=1;
        end
        function self = OptimizationObject(shape,name)
            % Create an MX symbolic object either by shape & name, or by
            % wrapping an existing MX object.
            create_new_object = true;
            if isa(shape, 'casadi.MX')
                assert(shape.is_symbolic(), 'Only an MX symbolic object can be wrapped.');
                the_object = shape;
                create_new_object = false;
            else
                if isscalar(shape)
                    shape = [shape 1];
                end
            end
            if create_new_object
                the_object = casadi.MX.sym(name,casadi.Sparsity.dense(shape(1),shape(2)));
            end
            self@casadi.MX(the_object);
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
           if size(self,1) == a && size(self,2)==b
              r = self;
           else
              r = casadi.reshape(self,a,b);
           end
        end
    end
    methods(Static, Hidden=true)
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

