classdef MyCallback < casadi.IterationCallback
    properties
      optisolve_instance
      mcallback
    end
        
    methods
        function self = MyCallback(optisolve_instance, mcallback)
          self.optisolve_instance = optisolve_instance;
          self.mcallback = mcallback;
        end
        function returncode = paren(self, solver)
            self.optisolve_instance.readoutputs(solver);
            self.mcallback(solver);
            returncode = 0;
        end
    end
end
