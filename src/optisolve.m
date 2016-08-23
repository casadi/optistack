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
            if length(varargin)>=1
                constraints = varargin{1};
            else
                constraints = {};
            end
            options = struct;
            if length(varargin)>=2
                options = varargin{2};
            end

            import casadi.*

            if ~iscell(constraints) || ~(isvector(constraints) || isempty(constraints))
                error('Constraints must be given as cell array: {x>=0,y<=0}');
            end
            if length(constraints)~=size(constraints,2)
                constraints = constraints';
            end
            

            [ gl_pure, gl_equality] = sort_constraints( constraints );
            [ scalar_objectives, twonorm_objectives, total_scalar_objective, total_objective ] = sort_objectives( objective );
            
            symbols = OptimizationObject.get_primitives({total_objective gl_pure{:}});
            
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
                g_helpers = {};
                for i = 1:length(gl_pure)
                   g_helpers = {g_helpers{:},MX.sym('g',gl_pure{i}.sparsity()) }; 
                end
                G_helpers = veccat(g_helpers{:});

                self.Ghelper = Function('Ghelper',{G_helpers},g_helpers);

                self.Ghelper_inv = Function('Ghelper_inv',g_helpers,{G_helpers});
            end
            
            codegen = false;
            if isfield(options,'codegen')
                codegen = options.codegen;
                options = rmfield(options,'codegen');
                options.jit = true;
            end
            
            opt = struct;
            
            gl_pure_v = MX();
            if ~isempty(gl_pure)
               gl_pure_v = veccat(gl_pure{:});
            end
            if ~isempty(twonorm_objectives)
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
            
            if isfield(options,'callback')
                mcallback = options.callback;
                options = rmfield(options,'callback');
                
                self.callback1 = MyCallback(self, mcallback);
                options.iteration_callback = self.callback1;
            end
            
            %opt.starcoloring_threshold = 1000;

            nlp = struct('x',X,'p',P,'f',total_objective,'g',gl_pure_v);

            self.solver = nlpsol('solver','ipopt', nlp, options);

            % Save to class properties
            self.symbols      = symbols;
            self.helper       = helper;
            self.helper_inv   = helper_inv;
            self.gl_equality  = gl_equality;
            
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
    
end
