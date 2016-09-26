x = optivar();
y = optivar();


p = optipar();
p.setValue(100);

nlp = optisolve((1-x)^2+p*(y-x^2)^2,{x^2+y^2<=1, x+y>=0});

optival(x)
optival(y)

assert(max(abs(optival(x)-0.7864151510041)<1e-9))
assert(max(abs(optival(y)-0.617698307517294)<1e-9))
assert(max(abs(optival(x^2+y^2)-1)<1e-7))

assert(max(abs(optival(p*y)-61.769830751729394)<1e-9))