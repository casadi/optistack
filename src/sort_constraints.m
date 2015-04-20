function [ gl_pure, gl_equality] = sort_constraints( gl )
    % Rewrites and determines nature of constraints, either g(x)<=0 or g(x)==0.
    % A user may write x>=y where x and y are variables.
    % In the `gl_pure` output, everything is brought to the left hand side
    % Parameters
    % ----------
    % gl : list of constraints, optional
    % Returns
    % -------
    % gl_pure : list of constraints in standard form
    % The constraints are rewritten as g(x)<=0 or g(x)==0
    % gl_equality : list of bools
    % For each entry in `gl_pure`, this list contains a boolean.

    gl_pure = {};
    gl_equality = [];
    for g = gl
        g = g{1};
        if g.isOperation(casadi.OP_LE) || g.isOperation(casadi.OP_LT)
        	gl_pure = {gl_pure{:},g.getDep(0) - g.getDep(1)};
            gl_equality = [gl_equality, false];
        elseif g.isOperation(casadi.OP_EQ)
        	gl_pure = {gl_pure{:},g.getDep(0) - g.getDep(1)};
            gl_equality = [gl_equality, true];
        else
            error('Constraint type unkown. Use ==, >= or <= .');
        end
    end
end

