classdef optisolve
    
    properties
            symbols
            Phelper_inv
            helper
            solver
            Ghelper_inv
            Ghelper
            helper_inv
            gl_equality
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
            import casadi.*

            if ~iscell(constraints)
                error('Constraints must be given as cell array: {x>=0,y<=0}');
            end

            [ gl_pure, gl_equality] = sort_constraints( constraints );
            symbols = OptimizationObject.get_primitives({objective gl_pure{:}});
            
            % helper functions for 'x'
            X = poor_veccat(symbols.x{:});
            helper = MXFunction({X},symbols.x);
            helper.init();

            helper_inv = MXFunction(symbols.x,{X});
            helper_inv.init();

            % helper functions for 'p' if applicable
            if isfield(symbols,'p')
              P = poor_veccat(symbols.p{:});

              self.Phelper_inv = MXFunction(symbols.p,{P});
              self.Phelper_inv.init();
              
            else
              P = [];
            end

            if ~isempty(gl_pure)
                g_helpers = {};
                for i = 1:length(gl_pure)
                   g_helpers = {g_helpers{:},MX.sym('g',gl_pure{i}.sparsity()) }; 
                end
                G_helpers = poor_veccat(g_helpers{:});

                self.Ghelper = MXFunction({G_helpers},g_helpers);
                self.Ghelper.init();

                self.Ghelper_inv = MXFunction(g_helpers,{G_helpers});
                self.Ghelper_inv.init();
            end

            nlp = MXFunction(nlpIn('x',X,'p',P), nlpOut('f',objective,'g',poor_veccat(gl_pure{:})));
            nlp.init();
            
            poor_veccat(gl_pure{:})

            self.solver = NlpSolver('ipopt', nlp);
            self.solver.init();

            % Save to class properties
            self.symbols      = symbols;
            self.helper       = helper;
            self.helper_inv   = helper_inv;
            self.gl_equality  = gl_equality;
            
            resolve(self);
        end
        
        function [] = resolve(self)
            % recall from class properties
            symbols      = self.symbols;
            helper       = self.helper;
            helper_inv   = self.helper_inv;
            gl_equality  = self.gl_equality;
            
            if ~isempty(gl_equality)
                % compose lbg
                for i=1:self.Ghelper_inv.getNumInputs()
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

            helper.setInput(self.solver.getOutput('x'));
            helper.evaluate();

            for i=1:length(symbols.x)
              v = symbols.x{i};
              v.setValue(full(helper.getOutput(i-1)));
            end

         end

    end
    
end
