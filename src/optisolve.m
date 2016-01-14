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
            helper = MXFunction('helper',{X},symbols.x);

            helper_inv = MXFunction('helper_inv',symbols.x,{X});

            % helper functions for 'p' if applicable
            if isfield(symbols,'p')
              P = veccat(symbols.p{:});

              self.Phelper_inv = MXFunction('Phelper_inv',symbols.p,{P});
              
            else
              P = MX.sym('p',0,1);
            end

            if ~isempty(gl_pure)
                g_helpers = {};
                for i = 1:length(gl_pure)
                   g_helpers = {g_helpers{:},MX.sym('g',gl_pure{i}.sparsity()) }; 
                end
                G_helpers = veccat(g_helpers{:});

                self.Ghelper = MXFunction('Ghelper',{G_helpers},g_helpers);

                self.Ghelper_inv = MXFunction('Ghelper_inv',g_helpers,{G_helpers});
            end
            
            codegen = false;
            if isfield(options,'codegen')
                codegen = options.codegen;
                options = rmfield(options,'codegen');
            end
            
            if isfield(options,'callback')
                mcallback = options.callback;
                options = rmfield(options,'callback');
                
                self.callback1 = MyCallback(self, mcallback);
                self.callback2 = Callback(self.callback1);
                options.iteration_callback = self.callback2;
            end
            
            opt = struct;
            if codegen
                opt.jit = true;
                opt.jit_options = struct('flags',char('-O3',''));
            end
            
            gl_pure_v = MX();
            if ~isempty(gl_pure)
               gl_pure_v = veccat(gl_pure{:});
            end
            if ~isempty(twonorm_objectives)
                F = [twonorm_objectives{:}];
                FF = MXFunction('nlp',{X,P},{F});

                JF = FF.jacobian();
                J_out = JF({X,P});
                J = J_out{1}';
                H = J*J';
                lam_f = MX.sym('lam_f');
                if isempty(scalar_objectives)
                    Hf = MXFunction('H',hessLagIn('x',X,'p',P,'lam_f',lam_f),hessLagOut('hess',lam_f*H),opt);
                else
                    lam_g = MX.sym('lam_g',size(gl_pure_v,1));
                    S = MXFunction('nlp',nlpIn('x',X,'p',P), nlpOut('f',total_scalar_objective,'g',gl_pure_v));
                    dS = S.derivative(0,1);
                    dS.printDimensions();
                    Hs = dS.jacobian(0,2,false,true);
                    Hs_out = Hs({X,P,lam_f,lam_g});
                    Hf = MXFunction('H',hessLagIn('x',X,'p',P,'lam_f',lam_f,'lam_g',lam_g),hessLagOut('hess',lam_f*H+Hs_out{1}),opt);
                end
                if isfield(options,'expand') && options.expand
                   Hf = SXFunction(Hf);
                end
                
                options.hess_lag = Hf;
            end
            
            nlp = MXFunction('nlp',nlpIn('x',X,'p',P), nlpOut('f',total_objective,'g',gl_pure_v),opt);

            if isfield(options,'expand') && options.expand
               nlp = nlp.expand();
            end
            self.solver = NlpSolver('solver','ipopt', nlp, options);

            % Save to class properties
            self.symbols      = symbols;
            self.helper       = helper;
            self.helper_inv   = helper_inv;
            self.gl_equality  = gl_equality;
            
            self.resolve();
        end
        
        function [] = jacspy(self)
            spy(sparse(DMatrix(self.solver.jacG().output(),1)))
        end
        
        function [] = resolve(self)
            % recall from class properties
            symbols      = self.symbols;
            helper       = self.helper;
            helper_inv   = self.helper_inv;
            gl_equality  = self.gl_equality;
            if ~isempty(gl_equality)
                % compose lbg
                for i=1:self.Ghelper_inv.nIn()
                    if gl_equality(i)
                        self.Ghelper_inv.setInput(0,i-1);
                    else
                        self.Ghelper_inv.setInput(-inf,i-1);
                    end
                end
                self.Ghelper_inv.evaluate();
                self.solver.setInput(self.Ghelper_inv.getOutput(),'lbg');
                self.solver.setInput(0,'ubg');
            end

            % compose lbx
            for i=1:length(symbols.x)
              helper_inv.setInput(symbols.x{i}.lb,i-1);
            end

            helper_inv.evaluate();
            self.solver.setInput(helper_inv.getOutput(),'lbx');  

            % compose x0
            for i=1:length(symbols.x)
              helper_inv.setInput(symbols.x{i}.init,i-1);
            end

            helper_inv.evaluate();
            self.solver.setInput(helper_inv.getOutput(),'x0'); 

            % compose ubx
            for i=1:length(symbols.x)
              helper_inv.setInput(symbols.x{i}.ub,i-1);
            end

            helper_inv.evaluate();
            self.solver.setInput(helper_inv.getOutput(),'ubx');    

            if isfield(symbols,'p')
                % compose p0
                for i=1:length(symbols.p)
                  self.Phelper_inv.setInput(symbols.p{i}.value,i-1);
                end

                self.Phelper_inv.evaluate();
                self.solver.setInput(self.Phelper_inv.getOutput(),'p');

            end
   
            self.solver.evaluate();
            self.readoutputs(self.solver);
            
         end
            
         function [] = readoutputs(self,solver)
            self.helper.setInput(solver.getOutput(0));
            self.helper.evaluate();

            for i=1:length(self.symbols.x)
              v = self.symbols.x{i};
              v.setValue(full(self.helper.getOutput(i-1)));
            end

         end

    end
    
end
