clear; clc; close all

% x = rand(12, 1);
% u = [zeros(3, 1); rand(3, 1)];
t = 0;
x = zeros(12, 1);
u = [zeros(3, 1); 1; 0; 0];

% [ts, xs] = ode45(@(t, x) dynamics(t, x, u), [0, 1], x);
x_dot = dynamics(t, x, u);
[A, B] = dynamics_jacobian(t, x, u);

function x_dot = dynamics(t, x, u)
[M,h,q_dot] = compute_dynamics(x(1:6), x(7:12));
x_dot = [q_dot; M\(u-h)];
end

function [A, B] = dynamics_jacobian(t, x, u)
[M,h,~] = compute_dynamics(x(1:6), x(7:12));
[dM,dh] = compute_dynamics_jacobian(x(1:6), x(7:12));
M_inv = inv(M);

A = zeros(6, 12);
for i=1:12
    A(:, i) = -M_inv*reshape(dM(:, i), [], 6)*M_inv*(u-h);
end
A = A+M\dh;

B = M_inv;
end