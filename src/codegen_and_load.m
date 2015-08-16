function out = codegen_and_load( name, fun)
import casadi.*

out = JitFunction('clang',fun);

end

