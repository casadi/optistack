N = 100;

% ode
ode = @(x,u) [u - x(1); x(1)];

% state: speed
S = optivar(N+1,1,'S');
% state: position
P = optivar(N+1,1,'P');

% control: 
U = optivar(N,1,'U');

% a variable to denote the final time
tf = optivar(1,1,'tf');
tf.setInit(1); % seconds

% Construct list of all constraints
g = {};

for k=1:N
   xk      = [S(k);   P(k)  ];
   xk_plus = [S(k+1); P(k+1)];
   
   % shooting constraint
   xf = rk4(ode,tf/N,xk,U(k));
   g = {g{:}, xk_plus==xf};
end


% path constraint
constr = @(P) 1-sin(2*pi*P)/2;

g = {g{:}, S <= constr(P)};

% Initialize with speed 1.
S.setInit(1);

U.setLb(0);
U.setUb(1);

g = {g{:}, S(1)==0, P(1)==0, P(end)==1};

optisolve(tf,g);

hold on
plot(optival(S));
plot(optival(P));
plot(constr(optival(P)),'r--');
stairs(optival(U),'b');
