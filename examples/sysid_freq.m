import casadi.*
close all

assert(false,'Requires some functionality from 2.4 which has not been transferred to 3.0 yet.')

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

param_truth = [5.625e-6;2.3e-4;1];
param_guess = [5;2;1];
scale = [1e-6;1e-4;1];

rng(1)

%%%%%%%%%%%% MODELING %%%%%%%%%%%%%%%%%%%%%
y  = MX.sym('y');
dy = MX.sym('dy');
u  = MX.sym('u');

states = [y;dy];
controls = u;

M = optivar();
c = optivar();
k = optivar();

params = [M;c;k];

rhs = [dy; (u-k*y-c*dy)/M];

% Form an ode function
ode = Function('ode',{states,controls,params},{rhs});

%%%%%%%%%%%% Creating a simulator %%%%%%%%%%
N_steps_per_sample = 10;
dt = 1/fs/N_steps_per_sample;

% Build an integrator for this system: Runge Kutta 4 integrator
k1 = ode(states,controls,params);
k2 = ode(states+dt/2.0*k1,controls,params);
k3 = ode(states+dt/2.0*k2,controls,params);
k4 = ode(states+dt*k3,controls,params);

states_final = states+dt/6.0*(k1+2*k2+2*k3+k4);

% Create a function that simulates one step propagation in a sample
one_step = Function('one_step',{states, controls, params},{states_final});

X = states;
for i=1:N_steps_per_sample
    X = one_step(X, controls, params);
end

% Create a function that simulates all step propagation on a sample
one_sample = Function('one_sample',{states, controls, params}, {X});

% speedup trick: expand into scalar operations
one_sample = one_sample.expand();

%%%%%%%%%%%% Simulating the system %%%%%%%%%%

all_samples = one_sample.mapaccum('all_samples', N);

% Choose an excitation signal
u_data = 0.1*rand(N,1);

x0 = DM([0,0]);
X_measured = all_samples(x0, u_data, repmat(param_truth,1,N));

y_data = X_measured(1,:)';

% Add some coloured noise
y_data = y_data + 0.2*filter(1,[1 -0.9],randn(N,1));

%%%%%%%%%%%% Identifying the simulated system: single shooting strategy %%%%%%%%%%

params_scale = params.*scale;

% Note, it is in general a good idea to scale your decision variables such
% that they are in the order of ~0.1..100
X_symbolic = all_samples(x0, u_data, repmat(params_scale,1,N));

e_time = y_data-X_symbolic(1,:)';

M.setInit(param_guess(1));
c.setInit(param_guess(2));
k.setInit(param_guess(3));

options = struct;
options.codegen = true;

disp('Single shooting without frequency information ...')

% Hand in a vector objective -> interpreted as 2-norm
% such t hat Gauss-Newton can be performed
optisolve(e_time,{},options);

disp('sol')
optival(M)*1e-6
optival(c)*1e-4
optival(k)

%rms(optival(e_time))
sqrt(mean(optival(e_time).^2)) % rms is part of the signal toolbox

%%

% Obtain the continuous time state space description symbolically
A = jacobian(rhs,states);
B = jacobian(rhs,controls);

C = jacobian(rhs(1),states);
D = jacobian(rhs(1),controls);

% And numerically
A_sol = optival(A);
B_sol = optival(B);
C_sol = optival(C);
D_sol = optival(D);

% Obtain those matrices for the true system as well
A_truth = optival(substitute(A,params,param_truth./scale));
B_truth = optival(substitute(B,params,param_truth./scale));
C_truth = optival(substitute(C,params,param_truth./scale));
D_truth = optival(substitute(D,params,param_truth./scale));

wplot = logspace(-2,1,1000);

% Continuous time model!
f_sol = squeeze(freqresp(ss(A_sol,B_sol,C_sol,D_sol),wplot,'rad/s'));
f_truth = squeeze(freqresp(ss(A_truth,B_truth,C_truth,D_truth),wplot,'rad/s'));

clf()
figure(1)
loglog(wplot,abs(f_truth))
hold on
loglog(wplot,abs(f_sol))


%% We want more weight on low frequencies

N_w = 100;
w_fit = logspace(-2,-1,N_w);

% Compute the 'true' frequency response
% In reality, we would obtain this through some fft
f_truth = squeeze(freqresp(ss(A_truth,B_truth,C_truth,D_truth),w_fit,'rad/s'));

% Construct a symbolic function that evaluates
% the magnitude of frequency reponse for one frequency
wsym = MX.sym('w');

% Note: solve should eventually be replaced by backslash (\),
% but this is not yet supported in casadi (as of 2.4.1)
Freal = @(R,I) D + C*solve(R + I*(solve(R,I,'lapacklu')),B,'lapacklu');
Fimag = @(R,I) -C*solve(I + R*(solve(I,R,'lapacklu')),B,'lapacklu');

Fmag = Function('Fmag',{wsym,params},{Freal(-A,wsym*eye(2)),Fimag(-A,wsym*eye(2))});

% Derive a function that evaluate Fmag for all frequencies
Fmag_all = Fmag.map('map',N_w);

[Freal, Fimag] = Fmag_all(w_fit,repmat(params_scale,1,N_w));

% The error between modeled frequency response and the truth
e_freq = [Freal' - real(f_truth);Fimag' - imag(f_truth)];

weight = 8;

e_total = [e_time;weight*e_freq];

options = struct;
% Note: codegen of linear solvers is missing in casadi still; that makes this slow
options.codegen = true;
options.tol = 1e-5;

optisolve(e_total,{},options);

disp('sol (with freq)')
optival(M)*1e-6
optival(c)*1e-4
optival(k)

A_solf = optival(A);
B_solf = optival(B);
C_solf = optival(C);
D_solf = optival(D);

f_sol = squeeze(freqresp(ss(A_solf,B_solf,C_solf,D_solf),wplot,'rad/s'));

loglog(wplot,abs(f_sol),'--')

legend('truth','fit (time)','fit (time+freq)')
