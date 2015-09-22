import casadi.*

% In this example, we fit a nonlinear model to measurements
%
% This example uses more advanced constructs than the vdp* examples:
% Since the number of control intervals is potentially very large here,
% we use memory-efficient Map and MapAccum, in combination with
% codegeneration.
%
% We will be working with a 2-norm objective:
% || y_measured - y_simulated ||_2^2
%
% This form is well-suited for the Gauss-Newton Hessian approximation.

%%%%%%%%%%% SETTINGS %%%%%%%%%%%%%%%%%%%%%
N = 10000;  % Number of samples
fs = 610.1; % Sampling frequency [hz]

param_truth = [5.625e-6;2.3e-4;1;4.69];
param_guess = [5;2;1;5];
scale = [1e-6;1e-4;1;1];

%%%%%%%%%%%% MODELING %%%%%%%%%%%%%%%%%%%%%
y  = MX.sym('y');
dy = MX.sym('dy');
u  = MX.sym('u');

states = [y;dy];
controls = u;

M = MX.sym('x');
c = MX.sym('c');
k = MX.sym('k');
k_NL = MX.sym('k_NL');

params = [M;c;k;k_NL];

rhs = [dy; (u-k_NL*y.^3-k*y-c*dy)/M];

% Form an ode function
ode = MXFunction('ode',{states,controls,params},{rhs});

%%%%%%%%%%%% Creating a simulator %%%%%%%%%%
N_steps_per_sample = 10;
dt = 1/fs/N_steps_per_sample;

% Build an integrator for this system: Runge Kutta 4 integrator
k1 = ode({states,controls,params});
k2 = ode({states+dt/2.0*k1{1},controls,params});
k3 = ode({states+dt/2.0*k2{1},controls,params});
k4 = ode({states+dt*k3{1},controls,params});

states_final = states+dt/6.0*(k1{1}+2*k2{1}+2*k3{1}+k4{1});

% Create a function that simulates one step propagation in a sample
one_step = MXFunction('one_step',{states, controls, params},{states_final});

X = states;
for i=1:N_steps_per_sample
    Xn = one_step({X, controls, params});
    X = Xn{1};
end

% Create a function that simulates all step propagation on a sample
one_sample = MXFunction('one_sample',{states, controls, params}, {X});

% speedup trick: expand into scalar operations
one_sample = one_sample.expand();

%%%%%%%%%%%% Simulating the system %%%%%%%%%%

all_samples = one_sample.mapaccum('all_samples', N);

% Choose an excitation signal
u_data = 0.1*rand(N,1);

x0 = DMatrix([0,0]);
all_samples_out = all_samples({x0, u_data, repmat(param_truth,1,N) });
X_measured = all_samples_out{1};

y_data = X_measured(1,:)';

%%%%%%%%%%%% Identifying the simulated system: single shooting strategy %%%%%%%%%%

% Note, it is in general a good idea to scale your decision variables such
% that they are in the order of ~0.1..100
all_samples_out = all_samples({x0, u_data, repmat(params.*scale,1,N) });
X_symbolic = all_samples_out{1};

e = y_data-X_symbolic(01,:)';

% just-in-time compilation for extra speedup
options = struct;
options.jit = true;
options.jit_options = struct('flags',char('-O3',''));

nlp = MXFunction('nlp', nlpIn('x',params), nlpOut('f',0.5*inner_prod(e,e)),options);

gradF = nlp.gradient();
jacG = nlp.jacobian('x','g');

gradF.derivative(0, 1);

J = jacobian(e,params);
sigma = MX.sym('sigma');
hessLag = MXFunction('H',hessLagIn('x',params,'lam_f',sigma),hessLagOut('hess',sigma*J'*J),options);

options = struct;
options.hess_lag = hessLag;
options.grad_f = gradF;
options.jac_g = jacG;
solver = NlpSolver('solver','ipopt', nlp, options);

sol = solver(struct('x0',param_guess));

sol.x.*scale

assert(max(full(abs(sol.x.*scale-param_truth)))<1e-8)
