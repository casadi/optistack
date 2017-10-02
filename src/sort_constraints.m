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
    gl_equality = false(1,0);
    for g = gl
        h = g{1};
        if h.is_op(casadi.OP_LE) || h.is_op(casadi.OP_LT)
            args = {};
            while h.is_op(casadi.OP_LE) || h.is_op(casadi.OP_LT)
               args = [args {h.dep(1)}]; %#ok<*AGROW>
               h = h.dep(0);
            end
            args = [args {h}];
            for i=1:length(args)-1
                gl_pure = [gl_pure, {args{i+1} - args{i}}];
                gl_equality = [gl_equality, false];
            end
        elseif h.is_op(casadi.OP_EQ)
        	gl_pure = [gl_pure, {h.dep(0) - h.dep(1)}];
            gl_equality = [gl_equality, true];
        else
            error('Constraint type unkown. Use ==, >= or <= .');
        end
    end
end

