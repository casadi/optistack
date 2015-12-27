function [ varargout ] = optival( varargin )
    import casadi.*
    symbols = OptimizationObject.get_primitives(varargin);
    
    hassymbols = false;
    symbolsx = {MX.sym('dummy')}; % bug in casadi typemaps: {} does not work
    if isfield(symbols,'x')
       symbolsx = symbols.x;
       hassymbols = true;
    end
    
    f = MXFunction('f',symbolsx,varargin);
    
    if hassymbols
        for i=1:length(symbolsx)
           666
           f.setInput(optival(symbolsx{i}),i-1); 
        end
    end
    
    f.evaluate();
    
    varargout = {};
    for i=1:length(varargin)
       varargout = {varargout{:},full(f.getOutput(i-1))};
    end
    
end

