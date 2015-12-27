n = 10;
m = 10000;

A = rand(m,n);
b = rand(m,1);

x = optivar(n,1,'x');

% This example shows a least-squares computation:
% min_x  || F(x) ||^2_2
%
% with F(x) = Ax+b
%

% A naive way to do this is:
tic
e = A*x+b;
optisolve(0.5*e'*e);
toc

x1 = optival(x);

% Solving with Gauss-Newton: pass vector as objective
tic
optisolve(A*x+b);
toc

x2 = optival(x);


% Solutions are the same
assert(max(abs(x1-x2))<1e-10)


