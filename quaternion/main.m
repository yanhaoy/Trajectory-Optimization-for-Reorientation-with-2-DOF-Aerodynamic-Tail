clear; clc; close all

% Initialize model
model = init_dynamics();

t = 0;
test=load('test');
% System states, quaternions from inertial frame of reference to body; 
% quaternions from body to tail
x = [1; zeros(3, 1); 1; zeros(3, 1); zeros(6, 1)];

% Input of tail's motor roll and pitch
u = [5; -1];

% Continuous dynamics
x_dot = dynamics(model, x, u);

% Continuous dynamics Jacobian
[A, B] = dynamics_jacobian(model, x, u);

% Discrete dynamics with explicit RK4
xk1 = discrete_dynamics(model, x, u, 0.02);

% Discrete dynamics Jacobian with explicit RK4
[Ak, Bk] = discrete_dynamics_jacobian(model, x, u, 1e-3);

% Simulation with Matlab ODE45
[ts, xs] = ode45(@(t, x) dynamics(model, x, u), [0, 1], x);

function x_dot = dynamics(model, x, u)

tmp = model.dynamics('x', x, 'u', u);
x_dot = full(tmp.x_dot);

end

function [A, B] = dynamics_jacobian(model, x, u)

tmp = model.dynamics_jacobian('x', x, 'u', u);
A = full(tmp.A);
B = full(tmp.B);

end

function xk1 = discrete_dynamics(model, x, u, dt)

tmp = model.discrete_dynamics('x', x, 'u', u, 'dt', dt);
xk1 = full(tmp.xk1);

end

function [Ak, Bk] = discrete_dynamics_jacobian(model, x, u, dt)

tmp = model.discrete_dynamics_jacobian('x', x, 'u', u, 'dt', dt);
Ak = full(tmp.Ak);
Bk = full(tmp.Bk);

end