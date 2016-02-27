classdef optisolve < handle
    
    properties
            symbols
            Phelper_inv
            helper
            solver
            Ghelper_inv
            Ghelper
            helper_inv
            gl_equality
            callback1
            callback2
            extra
            nx
            ng
            np
    end
    
    methods

        function [ self ] = optisolve( objective, varargin )
            %   optisolve(objective)
            %   optisolve(objective,constraints)
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
                    lam_g = MX.sym('lam_g',size(gl_pure_v,1));
                    S = Function('nlp',nlpIn('x',X,'p',P), nlpOut('f',total_scalar_objective));
                    dS = S.derivative(0,1);
                    Hs = dS.jacobian(0,2,false,true);
                    Hs_out = Hs({X,P,lam_f,0});
                    Hf = Function('nlp_hess_l',hessLagIn('x',X,'p',P,'lam_f',lam_f),hessLagOut('hess',lam_f*H+Hs_out{1}),opt);
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
            opt.monitor = char('eval_hess','eval_f','eval_grad_f','eval_jac_g');
            
            nlp = Function('nlp',struct('x',X,'p',P,'f',total_objective,'g',gl_pure_v),char('x','p'),char('f','g'),opt);

            if isfield(options,'expand') && options.expand
               nlp = nlp.expand();
               options = rmfield(options,'expand');
            end

            if codegen
                disp('Computing derivatives')
                grad_f = nlp.gradient();
                jac_g = nlp.jacobian(0,1);
                if isfield(options,'hess_lag')
                  hess_lag = options.hess_lag;
                else
                  grad_lag = nlp.derivative(0,1);
                  hess_lag = grad_lag.jacobian(0,2,false,true);
                  
                  hess_lag_ins = struct('x',X,'p',P,'lam_f',MX.sym('lam_f',hess_lag.sparsity_in(2)),'lam_g',MX.sym('lam_g',hess_lag.sparsity_in(3)));
                  hess_lag_ins2 = struct('der_x',X,'der_p',P);
                  hess_lag_ins2.adj0_f = hess_lag_ins.lam_f;
                  hess_lag_ins2.adj0_g = hess_lag_ins.lam_g;
                  
                  out = hess_lag.call(hess_lag_ins2);
                  hess_lag_ins.hess_gamma_x_x = triu(out.dadj0_x_dder_x);
                  hess_lag = Function('nlp_hess_l',hess_lag_ins,char('x','p','lam_f','lam_g'),char('hess_gamma_x_x'));
                  hess_lag.generate('nlp_hess_l');

                end
                disp('Codegenerating')
                nlp.generate('nlp');
                grad_f_ins = struct('x',X,'p',P);
                out = grad_f.call(grad_f_ins);
                grad_f_ins.grad_f_x = out.grad;
                grad_f_ins.f = out.f;
                grad_f = Function('nlp_grad_f',grad_f_ins,char('x','p'),char('f','grad_f_x'));
                grad_f.generate('nlp_grad_f');
                jac_g_ins = struct('x',X,'p',P);
                out = jac_g.call(jac_g_ins);
                jac_g_ins.jac_g_x = out.dg_dx;
                jac_g_ins.g = out.g;
                jac_g = Function('nlp_jac_g',jac_g_ins,char('x','p'),char('g','jac_g_x'));
                jac_g.generate('nlp_jac_g');
                hess_lag.generate('nlp_hess_l');
                
                jit_options = struct('flags',char('-O3',''));%,'plugin_libs',char('linearsolver_lapacklu',''));

                disp('Compiling')
                nlp_compiler = Compiler('nlp.c','clang',jit_options);
                nlp = external('nlp',nlp_compiler,struct);
                grad_f_compiler = Compiler('nlp_grad_f.c','clang',jit_options);
                grad_f = external('nlp_grad_f',grad_f_compiler,struct);
                jac_g_compiler = Compiler('nlp_jac_g.c','clang',jit_options);
                jac_g = external('nlp_jac_g',jac_g_compiler,struct);
                hess_lag_compiler = Compiler('nlp_hess_l.c','clang',jit_options);
                hess_lag = external('nlp_hess_l',hess_lag_compiler,struct);
                
                extra = struct;
                extra.nlp = nlp;
                extra.nlp_compiler = nlp_compiler;
                extra.grad_f = grad_f;
                extra.grad_f_compiler = grad_f_compiler;
                extra.jac_g = jac_g;
                extra.jac_g_compiler = jac_g_compiler;
                extra.hess_lag = hess_lag;
                extra.hess_lag_compiler = hess_lag_compiler;
                
                self.extra = extra;
                
                options.grad_f = grad_f;
                options.jac_g = jac_g;
                options.hess_lag = hess_lag;
             
            end


            self.solver = nlpsol('solver','ipopt', nlp, options);

            % Save to class properties
            self.symbols      = symbols;
            self.helper       = helper;
            self.helper_inv   = helper_inv;
            self.gl_equality  = gl_equality;
            
            self.resolve();
        end
        
        function [] = jacspy(self)
            spy(sparse(DM(self.solver.jacG().output(),1)))
        end
        
        function [] = resolve(self)
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
            
         function [] = readoutputs(self,solver_out)
            helper_outputs = self.helper.call({solver_out.x});

            for i=1:length(self.symbols.x)
              v = self.symbols.x{i};
              v.setValue(full( helper_outputs{i}));
            end

         end

    end
    
end
