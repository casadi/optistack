assert(false,'Requires some funcitonality from 2.4 which has not been transferred to 3.0 yet.')

global x
global y

x = optivar();
y = optivar();

hold on
xlim([-1,1])
ylim([-1,1])
opts.callback = @mycallback;

nlp = optisolve((1-x)^2+100*(y-x^2)^2,{x^2+y^2<=1, x+y>=0},opts);

optival(x)
optival(y)

assert(max(abs(optival(x)-0.7864151510041)<1e-9))
assert(max(abs(optival(y)-0.617698307517294)<1e-9))
assert(max(abs(optival(x^2+y^2)-1)<1e-7))
