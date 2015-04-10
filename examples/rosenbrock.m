x = optivar();
y = optivar();

optisolve((1-x)^2+100*(y-x^2)^2,[x^2+y^2<=1, x+y>=0]);

value(x)
value(y)
optisolve