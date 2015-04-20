x = optivar();
y = optivar();

optisolve((1-x)^2+100*(y-x^2)^2,{x^2+y^2<=1, x+y>=0})

value(x)
value(y)

assert(max(abs(value(x)-0.7864151510041)<1e-9))
assert(max(abs(value(y)-0.617698307517294)<1e-9))