A = optivar(3,3);

% inverse
B = [1 0 -2; 3 1 -1;2 1 0];

optisolve(sum_square(A*B-eye(3)));

As = optival(A)

assert(max(max(abs(As-inv(B)))<1e-7));

A = optivar(4,3);
% moore-penrose pseudo-inverse
B = [1 0 -2 3; 3 1 -1 1;2 1 0 9];

optisolve(sum_square(A*B-eye(4)));

As = optival(A)

assert(max(max(abs(As-pinv(B)))<1e-7));






