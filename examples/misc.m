x = optivar();
y = optivar();

optisolve((x-1)^2+x);

optival(x)

assert(max(abs(optival(x)-0.5)<1e-9))

%%

p = optipar();

p.setValue(100);

nlp = optisolve((1-x)^2+p*(y-x^2)^2,{x^2+y^2<=1, x+y>=0});

optival(x)
optival(y)

assert(max(abs(optival(x)-0.7864151510041)<1e-9))
assert(max(abs(optival(y)-0.617698307517294)<1e-9))
assert(max(abs(optival(x^2+y^2)-1)<1e-7))

p.setValue(10);
nlp.resolve();

assert(max(abs(optival(x)-0.788740490289645)<1e-9))
assert(max(abs(optival(y)-0.614726302810607)<1e-9))
assert(max(abs(optival(x^2+y^2)-1)<1e-7))

%%
nlp = optisolve((1-x)^2+100*(y-x^2)^2,{x^2+y^2==1});

optival(x)
optival(y)

assert(max(abs(optival(x)-0.7864151510041)<1e-8))
assert(max(abs(optival(y)-0.617698307517294)<1e-8))
assert(max(abs(optival(x^2+y^2)-1)<1e-7))
