close all
clear variables
clc

import casadi.*

N = 30;  % Control discretization
T = 3.0; % End time
nx = 1;
nu = 1;

% Declare optimization variables
U  = optivar(nu,N);    % control
X  = optivar(nx,N+1);  % states

% System dynamics
ode = @(x,u) [(1 - x(1))*x(1)  + u(1)];

% Collection of constraints
g = {};

% Objective function
J=0;

x0 = optipar(1,1,'x0');
for i = 1:N
    % Propagate state
    x_next = rk4(ode,T/N,X(:,i),U(:,i));
    
    % Gap closing constraint
    g = {g{:}, X(:,i+1) == x_next};
    
    % Construct an objective
    J = J + x_next(1)^2 + U(1,i)^2;
end

% Terminal constraint
g = {g{:} X(1,end) == 0 };
% Initial constraint
g = {g{:} X(1,1) == x0 };

% Path constraints
g = {g{:} -1 <= X <= 1 };

% Initialization
X.setInit(0);

x0.setValue(0.01);
sol = optisolve(J,g);

x = 4;

N_mpc = 50;

Xall = [];
Uall = [];

figure()

while true
    x0.setValue(x);
    sol.resolve();
    
    Usol = optival(U);
    u = Usol(:,1);
    
    Xall = [Xall x];
    Uall = [Uall u];
    
    % Do not remember history larger than 50
    if size(Xall,2)>50
       Xall = Xall(:,2:50);
       Uall = Uall(:,2:50);
    end
    
    % Bring the physical system one step further
    x = rk4(ode,T/N,x,u);
    
    % Perturb if ctrl is pressed
    modifiers = get(gcf,'currentModifier');
    if ismember('control',modifiers)
       x=x+1; 
    end
    
    clf
    title('Press <ctrl> to perturb')
    hold on
    plot(Xall,'b','LineWidth',3);
    stairs(Uall,'r','LineWidth',3);
    
    plot(size(Xall,2):size(Xall,2)+N,optival(X),'b');
    stairs(size(Xall,2):size(Xall,2)+N-1,optival(U),'r');
    
    grid on
    
    pause(0.01);
end

