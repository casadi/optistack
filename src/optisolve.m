classdef optisolve < handle
    % OPTISOLVE  Solve an OPTISTACK NLP
    % The optisolve class performs the transcription of
    % an expression graph to an NLP, and solves it by Ipopt.
    %
    % Example:
    %  x = optivar();
    %  y = optivar();
    %  nlp = optisolve((1-x)^2+100*(y-x^2)^2,{x^2+y^2<=1, x+y>=0});
    %  optival(x)
    %  optival(y)
    %
    % optisolve Methods:
    %    optisolve - Solve an NLP
    %    resolve - Solve an NLP again; use together with an updated optipar

    properties(Access=?MyCallback)
            symbols
            Phelper_inv
            helper
            Ghelper_inv
            Ghelper
            helper_inv
            gl_equality
            callback1
            callback2
            extra
    end

    properties(SetAccess=private)
            solver % CasADi solver object
            nx % Number of decision variable
            ng % Number of constraints
            np % Number of parameters
    end
    
    % Properties used for fmincon
    properties(SetAccess=private)
        m_objfunc   % The objective function, as a Function object
        m_objfunc_g % The gradient of objective function, that returns both f and grad_f
        m_ieqc      % Inequality constraints g <= 0
        m_ieqc_j    % Jacobian of inequality constraints (and the constraints)
        m_eqc       % Equality constraints g == 0, including jacobian
        m_eqc_j     % Jacobian of equality constraints (and the constraints)
        
        m_p0        % The values of parameters in the current resolve()
        
        m_fmincon_opt   % Options for fmincon()
        
        m_fmincon = false;  % If fmincon is used
        
        fmincon_flag = NaN; % exit flag of fmincon, after it has solved the problem
    end

    methods

        function [ self ] = optisolve( objective, varargin )
            % OPTISOLVE Solve an OPTISTACK NLP
            %   optisolve(objective)
            %   optisolve(objective,constraints)
            %   optisolve(objective,constraints,options)
            %
            %     objective   - expression graph;
            %                   constructed from performing mathematical
            %                   operations on optivar and optipar instances
            %                   Examples:
            %                       sin(x)
            %                       x'*x
            %
            %     constraints - cell array of expression graphs
            %                   Examples:
            %                        {x>=0}
            %                        {x+y==0,0<=x<=1}
            %
            %     options     - struct
            %                   common options:
            %                    expand (true|false) trade speed for memory
            %                    codegen (true|false) compile-and-load code
            %                    ipopt (struct) pass options on to ipopt
            %
            %     solver      - string
            %                   the solver name, default: ipopt
            %
            % The solution of the NLP may be obtained using optival:
            %
            %  x = optivar();
            %  y = optivar();
            %  nlp = optisolve((1-x)^2+100*(y-x^2)^2,{x^2+y^2<=1, x+y>=0});
            %  optival(x)
            %  optival(y)
            %
            %
            % See also OPTIVAR, OPTIPAR, OPTIVAL
            
            % optisolve_start = tic;
            
            if length(varargin) >= 1
                constraints = varargin{1};
            else
                constraints = {};
            end
            options = struct;
            if length(varargin) >= 2 && isstruct(varargin{2})
                options = varargin{2};
            end
            solver = 'ipopt';
            if length(varargin) >= 3
                assert(ischar(varargin{3}), 'Solver name must be a string.');
                solver = varargin{3};
            end
            self.m_fmincon = strcmpi(solver, 'fmincon');
            
            import casadi.*
            
            if ~iscell(constraints) || ~(isvector(constraints) || isempty(constraints))
                error('Constraints must be given as cell array: {x>=0,y<=0}');
            end
            if length(constraints)~=size(constraints,2)
                constraints = constraints';
            end
            
            
            [ gl_pure, gl_equality] = sort_constraints( constraints );
            [ scalar_objectives, twonorm_objectives, total_scalar_objective, total_objective ] = sort_objectives( objective );
            
            symbols = OptimizationObject.get_primitives([{total_objective} gl_pure]);
            
            % helper functions for 'x'
            X = veccat(symbols.x{:});
            helper = Function('helper',{X},symbols.x);
            
            helper_inv = Function('helper_inv',symbols.x,{X});
            
            % helper functions for 'p' if applicable
            if isfield(symbols,'p')
                P = veccat(symbols.p{:});
                
                self.Phelper_inv = Function('Phelper_inv',symbols.p,{P});
                
            else
                P = MX.sym('p',0,1);
            end
            
            self.np = size(P,1);
            self.nx = size(X,1);
            
            if ~isempty(gl_pure)
                g_helpers = cell(1, length(gl_pure));
                for k = 1:length(gl_pure)
                    g_helpers{k} = MX.sym('g',gl_pure{k}.sparsity());
                end
                G_helpers = veccat(g_helpers{:});
                
                self.Ghelper = Function('Ghelper',{G_helpers},g_helpers);
                
                self.Ghelper_inv = Function('Ghelper_inv',g_helpers,{G_helpers});
            end
            
            opt = struct;
            
            gl_pure_v = MX();
            if ~isempty(gl_pure)
                gl_pure_v = veccat(gl_pure{:});
            end
            if ~isempty(twonorm_objectives) && ~self.m_fmincon
                F = [twonorm_objectives{:}];
                FF = Function('nlp',{X,P},{F});
                
                JF = FF.jacobian();
                J_out = JF.call({X,P});
                J = J_out{1}';
                H = J*J';
                lam_f = MX.sym('lam_f');
                if isempty(scalar_objectives)
                    Hf = Function('nlp_hess_l',struct('x',X,'p',P,'lam_f',lam_f,'hess_gamma_x_x',lam_f*triu(H)),char('x', 'p', 'lam_f', 'lam_g'),char('hess_gamma_x_x'),opt);
                else
                    S = Function('nlp',struct('x',X,'p',P,'f',total_scalar_objective),char('x','p'),char('f','g'));
                    dS = S.derivative(0,1);
                    Hs = dS.jacobian(0,2,false,true);
                    Hs_out = Hs.call({X,P,lam_f,0});
                    Hf = Function('nlp_hess_l',struct('x',X,'p',P,'lam_f',lam_f,'hess_gamma_x_x',triu(lam_f*H+Hs_out{1})),char('x', 'p', 'lam_f', 'lam_g'),char('hess_gamma_x_x'),opt);
                end
                if isfield(options,'expand') && options.expand
                    Hf = Hf.expand();
                end
                
                options.hess_lag = Hf;
            end
            self.ng = size(gl_pure_v,1);
            
            %opt.starcoloring_threshold = 1000;
            
            if self.m_fmincon
                % If the solver is fmincon, save the objective function and
                % constraints.
                
                [self.m_fmincon_opt, ~] = self.fmincon_process_options(options);
                
                self.m_objfunc = Function('obj', {X, P}, {total_objective}, {'x', 'p'}, {'f'}, opt);
                self.m_objfunc_g = self.m_objfunc.gradient('x', 'f');
                
                % In the following, we must use find() because CasADi does not
                % support logical indexing (at the moment).
                if any(gl_equality)
                    % If there is any equality constraint
                    self.m_eqc = Function('eqc', {X, P}, {gl_pure_v(find(gl_equality))}, {'x', 'p'}, {'g'}, opt); %#ok<FNDSB>
                    self.m_eqc_j = self.m_eqc.jacobian('x', 'g');
                else
                    self.m_eqc = [];
                    self.m_eqc_j = [];
                end
                
                if ~all(gl_equality)
                    % If there is any inequality constraint
                    self.m_ieqc = Function('ieqc', {X, P}, {gl_pure_v(find(~gl_equality))}, {'x', 'p'}, {'g'}, opt); %#ok<FNDSB>
                    self.m_ieqc_j = self.m_ieqc.jacobian('x', 'g');
                else
                    self.m_ieqc = [];
                    self.m_ieqc_j = [];
                end
                
            else
                % Use CasADi's supported solver
                if isfield(options,'codegen')
                    codegen = options.codegen;
                    options = rmfield(options,'codegen');
                    options.jit = codegen;
                end
                
                if isfield(options,'quiet')
                    quiet = options.quiet;
                    options = rmfield(options,'quiet');
                    if quiet
                        options.print_time = false;
                        switch solver
                            case 'ipopt'
                                options.ipopt.print_level = 0;
                            otherwise
                        end
                    end
                end
                
                if isfield(options,'callback')
                    mcallback = options.callback;
                    options = rmfield(options,'callback');
                    
                    self.callback1 = MyCallback(self, mcallback);
                    options.iteration_callback = self.callback1;
                end
                
                nlp = struct('x',X,'p',P,'f',total_objective,'g',gl_pure_v);
                
                self.solver = nlpsol('solver',solver, nlp, options);
            end
            
            % Save to class properties
            self.symbols      = symbols;
            self.helper       = helper;
            self.helper_inv   = helper_inv;
            self.gl_equality  = gl_equality;
            
            % fprintf('optisolve overhead = %g sec;', toc(optisolve_start));
            
            self.resolve();
        end

        function [] = resolve(self)
            % Solve NLP again; use together with an updated optipar
            %
            % Example:
            % x=optivar(3,1);
            %
            % r=optipar();
            % r.setValue(1)
            % sol = optisolve([1 1 0]*x,{x'*x==r});
            % optival(x) % [-0.7071;-0.7071;0]
            %
            % r.setValue(2)
            % sol.resolve();
            %
            % optival(x) % [-1;-1;0]
            
            if self.m_fmincon, self.resolve_fmincon(); return; end

            % recall from class properties
            symbols      = self.symbols;
            helper       = self.helper;
            helper_inv   = self.helper_inv;
            gl_equality  = self.gl_equality;
            solver_inputs = struct;
            if ~isempty(gl_equality)
                Ghelper_inv_inputs = cell(1,self.Ghelper_inv.n_in);
                % compose lbg
                for i=1:self.Ghelper_inv.n_in()
                    if gl_equality(i)
                        Ghelper_inv_inputs{i} = 0;
                    else
                        Ghelper_inv_inputs{i} = -inf;
                    end
                end
                out = self.Ghelper_inv.call(Ghelper_inv_inputs);
                solver_inputs.lbg = out{1};
                solver_inputs.ubg = 0;
            end

            helper_inv_inputs = cell(1,helper_inv.n_in);
            % compose lbx
            for i=1:length(symbols.x)
              helper_inv_inputs{i} = symbols.x{i}.lb;
            end

            out = helper_inv.call(helper_inv_inputs);
            solver_inputs.lbx = out{1};

            helper_inv_inputs = cell(1,helper_inv.n_in);
            % compose x0
            for i=1:length(symbols.x)
              helper_inv_inputs{i} = symbols.x{i}.init;
            end

            out = helper_inv.call(helper_inv_inputs);
            solver_inputs.x0 = out{1};

            helper_inv_inputs = cell(1,helper_inv.n_in);
            % compose ubx
            for i=1:length(symbols.x)
              helper_inv_inputs{i} = symbols.x{i}.ub;
            end

            out = helper_inv.call(helper_inv_inputs);
            solver_inputs.ubx = out{1};

            if isfield(symbols,'p')
                Phelper_inv_inputs = cell(1,self.Phelper_inv.n_in);

                % compose p0
                for i=1:length(symbols.p)
                  Phelper_inv_inputs{i} = symbols.p{i}.value;
                end
                out = self.Phelper_inv.call(Phelper_inv_inputs);
                solver_inputs.p = out{1};
            end
            out = self.solver.call(solver_inputs);
            self.readoutputs(out);

         end
    end
    methods(Access=?MyCallback)
         function [] = readoutputs(self,solver_out)
            helper_outputs = self.helper.call({solver_out.x});

            for i=1:length(self.symbols.x)
              v = self.symbols.x{i};
              v.setValue(full( helper_outputs{i}));
            end
            
         end
         
    end
    
    methods(Access=protected)
        function [f, g] = fmincon_objfun(self, x)
            % Objective function for fmincon
            
            % Call the saved objective function and its gradient
            if nargout > 1
                out = self.m_objfunc_g.call(struct('x', x, 'p', self.m_p0));
                f = full(out.f);
                g = full(out.grad);
            else
                f = full(self.m_objfunc(x, self.m_p0));
            end
        end
        
        function [c, ceq, gradc, gradceq] = fmincon_nonlcon(self, x)
            % Constraint function for fmincon
            
            % Call the saved constraint functions and their gradients
            if nargout > 2
                inputs = struct('x', x, 'p', self.m_p0);    % inputs to function calls
            
                if isempty(self.m_ieqc_j)
                    c = [];
                    gradc = [];
                else
                    out_c = self.m_ieqc_j.call(inputs);
                    c = full(out_c.g);
                    if nargout > 2
                        gradc = full(out_c.dg_dx).';
                    end
                end
                
                if isempty(self.m_eqc_j)
                    ceq = [];
                    gradceq = [];
                else
                    out_ceq = self.m_eqc_j.call(inputs);
                    ceq = full(out_ceq.g);
                    if nargout > 2
                        gradceq = full(out_ceq.dg_dx).';
                    end
                end
            else
                if isempty(self.m_ieqc)
                    c = [];
                else
                    c = full(self.m_ieqc(x, self.m_p0));
                end
                
                if isempty(self.m_eqc)
                    ceq = [];
                else
                    ceq = full(self.m_eqc(x, self.m_p0));
                end
            end
        end
        
        function resolve_fmincon(self)
            % Solve NLP using fmincon
            % The functions involved in the optimization must have already
            % been saved by the constructor (the solver must be set to
            % 'fmincon' in the constructor).
            
            % recall from class properties
            symbols      = self.symbols;
            helper       = self.helper;
            helper_inv   = self.helper_inv;
            gl_equality  = self.gl_equality;
            
            solver_inputs = struct;

            helper_inv_inputs = cell(1,helper_inv.n_in);
            % compose lbx
            for i=1:length(symbols.x)
                helper_inv_inputs{i} = symbols.x{i}.lb;
            end
            
            out = helper_inv.call(helper_inv_inputs);
            solver_inputs.lbx = full(out{1});
            
            helper_inv_inputs = cell(1,helper_inv.n_in);
            % compose x0
            for i=1:length(symbols.x)
                helper_inv_inputs{i} = symbols.x{i}.init;
            end
            
            out = helper_inv.call(helper_inv_inputs);
            solver_inputs.x0 = full(out{1});
            
            helper_inv_inputs = cell(1,helper_inv.n_in);
            % compose ubx
            for i=1:length(symbols.x)
                helper_inv_inputs{i} = symbols.x{i}.ub;
            end
            
            out = helper_inv.call(helper_inv_inputs);
            solver_inputs.ubx = full(out{1});
            
            if isfield(symbols,'p')
                Phelper_inv_inputs = cell(1,self.Phelper_inv.n_in);
                
                % compose p0
                for i=1:length(symbols.p)
                    Phelper_inv_inputs{i} = symbols.p{i}.value;
                end
                out = self.Phelper_inv.call(Phelper_inv_inputs);
                self.m_p0 = out{1};
            end
            
            % disp(self.m_fmincon_opt);
            
            % Call fmincon()
            solver_outputs = struct();
            [solver_outputs.x, ~, self.fmincon_flag] = fmincon(@self.fmincon_objfun, solver_inputs.x0,...
                [], [], [], [], ...
                solver_inputs.lbx, solver_inputs.ubx, ...
                @self.fmincon_nonlcon, ...
                self.m_fmincon_opt);
            
            self.readoutputs(solver_outputs);
            
        end
    end
    
    methods(Static, Hidden=true)
        function [opt, options] = fmincon_process_options(options)
            % Process the options for fmincon(), returns remaining unused
            % options.
            if isfield(options, 'fmincon')
                % Save the options for fmincon
                opt = options.fmincon;
                options = rmfield(options, 'fmincon');
            else
                opt = optimoptions('fmincon');
            end
            
            if isfield(options,'quiet')
                quiet = options.quiet;
                options = rmfield(options,'quiet');
                if quiet
                    opt = optimoptions(opt, 'Display', 'off');
                end
            end
            
            % Set the options for gradient computations
            opt = optimoptions(opt, 'GradObj', 'on', 'GradConstr', 'on');
        end
    end
end
