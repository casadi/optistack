function [ varargout ] = optival( varargin )
    import casadi.*
    symbols = OptimizationObject.get_primitives(varargin);
    
    hassymbols = false;
    symbolsx = {MX.sym('dummy')}; % bug in casadi typemaps: {} does not work
    if isfield(symbols,'x')
       symbolsx = symbols.x;
       hassymbols = true;
    end
    
    f = Function('f',symbolsx,varargin);
    
    f_inputs = cell(1,f.n_in);
    if hassymbols
        for i=1:length(symbolsx)
           f_inputs{i} = optival(symbolsx{i}); 
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

