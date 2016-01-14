x = optivar();
y = optivar();

nlp = optisolve({(1-x)^2+100*(y-x^2)^2,[x;y]},{});

num1 = optival([x;y]);


e = [x;y];

f = 0.5*inner_prod(e,e);

nlp = optisolve((1-x)^2+100*(y-x^2)^2+f,{});

num2 = optival([x;y]);

assert(max(abs(num1-num2))<1e-9);

