function [ varargout ] = optival( varargin )
    import casadi.*
    symbols = OptimizationObject.get_primitives(varargin);
    
    f = MXFunction('f',symbols.x,varargin);
    f.init();
    
    for i=1:length(symbols.x)
       f.setInput(optival(symbols.x{i}),i-1); 
    end
    
    f.evaluate();
    
    varargout = {};
    for i=1:length(varargin)
       varargout = {varargout{:},full(f.getOutput(i-1))};
    end
    
end

