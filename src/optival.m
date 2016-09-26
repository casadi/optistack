function [ varargout ] = optival( varargin )
    % Evaluate an expression at the optimimum found by OPTISTACK.
    %
    % x_opt = optival(x);
    % e_opt = optival(x'*x);
    
    import casadi.*
    symbols_struct = OptimizationObject.get_primitives(varargin);
    
    symbols = {};
    hassymbols = false;
    if isfield(symbols_struct,'x')
       symbols = symbols_struct.x;
       hassymbols = true;
    end
    if isfield(symbols_struct,'p')
       symbols = [symbols symbols_struct.p];
       hassymbols = true;
    end

    if ~hassymbols
       symbols = {MX.sym('dummy')}; % bug in casadi typemaps: {} does not work
    end

    
    f = Function('f',symbols,varargin);
    
    f_inputs = cell(1,f.n_in);
    if hassymbols
        for i=1:length(symbols)
           f_inputs{i} = optival(symbols{i}); 
        end
    else
        f_inputs{1} = 0;
    end
    
    f_out = f.call(f_inputs);
    varargout = {};
    for i=1:length(varargin)
       varargout = {varargout{:},full(f_out{i})};
    end
    
end

