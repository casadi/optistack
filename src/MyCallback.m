classdef MyCallback < casadi.Callback
    properties
      optisolve_instance
      mcallback
    end
        
    methods
        function self = MyCallback(optisolve_instance, mcallback)
          self@casadi.Callback();
          self.optisolve_instance = optisolve_instance;
          self.mcallback = mcallback;
          opts = struct;
          opts.output_scheme = char('ret');
          opts.input_scheme = casadi.nlpsol_out();
          construct(self,'MyCallback', opts);
        end
        function returncode = eval(self, solver)
            solver_struct = struct;
            for i=1:casadi.nlpsol_n_out()
                solver_struct.(casadi.nlpsol_out(i-1)) = solver{i};
            end
            self.optisolve_instance.readoutputs(solver_struct);
            self.mcallback(solver_struct);
            returncode = {0};
        end
        function out = get_input_shape(self,i)
          n = casadi.nlpsol_out(i);
          if strcmp(n,'f')
            out = [1 1];
          elseif strcmp(n,'lam_x') || strcmp(n,'x')
            out = [self.optisolve_instance.nx 1];
          elseif strcmp(n,'lam_g') || strcmp(n,'g')
            out = [self.optisolve_instance.ng 1];
          elseif strcmp(n,'p')  || strcmp(n,'lam_p')
            out = [self.optisolve_instance.np 1];
          else
            out = [0 0];
          end
        end
        function out = get_n_in(self)
          out = casadi.nlpsol_n_out();
        end
        function out = get_n_out(self)
          out = 1;
        end
        
    end
end
