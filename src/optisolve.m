function [ optimal ] = optisolve( objective, constraints )
    import casadi.*
    
    if ~iscell(constraints)
        error('Constraints must be given as cell array: {x>=0,y<=0}');
    end

    [ gl_pure, gl_equality] = sort_constraints( constraints );
    symbols = OptimizationObject.get_primitives({objective gl_pure{:}});
    
    helper = MXFunction({vertcat(symbols.x{:})},symbols.x);
    helper.init();
    
    nlp = MXFunction(nlpIn('x',vertcat(symbols.x{:})), nlpOut('f',objective,'g',vertcat(gl_pure{:})));
    nlp.init();

    solver = NlpSolver('ipopt', nlp);
    solver.init();
    
    ng = length(gl_pure);
    
    lbg = zeros(ng,1);
    ubg = zeros(ng,1);

    for i = 1:length(gl_equality)
        if gl_equality(i)
            lbg(i) = 0;
            ubg(i) = 0;
        else
            lbg(i) = -inf;
            ubg(i) = 0;
        end

    solver.setInput(lbg,'lbg');
    solver.setInput(ubg,'ubg');
    solver.evaluate();

    helper.setInput(solver.getOutput('x'));
    helper.evaluate();
    
    for i=1:length(symbols.x)
      v = symbols.x{i};
      v.setValue(full(helper.getOutput(i-1)));
    end

end

