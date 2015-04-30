function [ optimal ] = optisolve( objective, constraints )
    import casadi.*
    
    if ~iscell(constraints)
        error('Constraints must be given as cell array: {x>=0,y<=0}');
    end

    [ gl_pure, gl_equality] = sort_constraints( constraints );
    symbols = OptimizationObject.get_primitives({objective gl_pure{:}});
    
    X = vertcat(symbols.x{:});
    helper = MXFunction({X},symbols.x);
    helper.init();
    
    helper_inv = MXFunction(symbols.x,{X});
    helper_inv.init();
    
    
    ng = length(gl_pure);
    
    g_helpers = {};
    for i = 1:ng
       g_helpers = {g_helpers{:},MX.sym('g',gl_pure{i}.sparsity()) }; 
    end
    G_helpers = vertcat(g_helpers{:});
    
    Ghelper = MXFunction({G_helpers},g_helpers);
    Ghelper.init();
    
    Ghelper_inv = MXFunction(g_helpers,{G_helpers});
    Ghelper_inv.init();
    
    
    nlp = MXFunction(nlpIn('x',X), nlpOut('f',objective,'g',vertcat(gl_pure{:})));
    nlp.init();

    solver = NlpSolver('ipopt', nlp);
    solver.init();

    nx = length(X);
    
    lbg = zeros(ng,1);
    ubg = zeros(ng,1);

        % compose lbg
    for i=1:length(gl_pure)
        if gl_equality(i)
            Ghelper_inv.setInput(0,i-1);
        else
            Ghelper_inv.setInput(-inf,i-1);
        end
    end
    
    Ghelper_inv.evaluate();
    solver.setInput(Ghelper_inv.getOutput(),'lbg');    
    
    % compose lbx
    for i=1:length(symbols.x)
      helper_inv.setInput(symbols.x{i}.lb,i-1);
    end
    
    helper_inv.evaluate();
    solver.setInput(helper_inv.getOutput(),'lbx');    
    
    % compose ubx
    for i=1:length(symbols.x)
      helper_inv.setInput(symbols.x{i}.ub,i-1);
    end
    
    helper_inv.evaluate();
    solver.setInput(helper_inv.getOutput(),'ubx');    
    

    solver.setInput(0,'ubg');
    solver.evaluate();

    helper.setInput(solver.getOutput('x'));
    helper.evaluate();
    
    for i=1:length(symbols.x)
      v = symbols.x{i};
      v.setValue(full(helper.getOutput(i-1)));
    end

end

