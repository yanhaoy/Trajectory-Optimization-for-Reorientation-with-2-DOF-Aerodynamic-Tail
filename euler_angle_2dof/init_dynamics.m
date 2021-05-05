function [model] = init_dynamics()

% addpath('..\utils\')
% addpath('..\utils\casadi\')
addpath('utils')
addpath('utils\casadi')
import casadi.*

% Body mass
m_b = 1;

% Tail mass
m_t = 1;

% Body inertia
J_b = [1; 1; 1];

% Tail inertia
J_t = [1; 1; 1];

% Translation from body to joint in body frame
l_b = [1; 0; 0];

% Translation from tail to joint in tail frame
l_t = [- 1; 0; 0];

% Generalized coordinates: Euler angles from inertial frame of reference to
% body, roll, pitch, and yaw; Euler angles from body to tail, roll and
% pitch, we don't use yaw here to keep consistency to the 2 motors
q = SX.sym('q', [5, 1]);

% Euler angle rates
q_dot = SX.sym('q_dot', [5, 1]);

%% Rotation matrices

% Inertial frame of reference to body frame
R_Ob = [cos(q(2))*cos(q(3)), -cos(q(2))*sin(q(3)), sin(q(2));
    cos(q(1))*sin(q(3)) + cos(q(3))*sin(q(1))*sin(q(2)), cos(q(1))*cos(q(3)) - sin(q(1))*sin(q(2))*sin(q(3)), -cos(q(2))*sin(q(1));
    sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(3))*sin(q(2)), cos(q(3))*sin(q(1)) + cos(q(1))*sin(q(2))*sin(q(3)),  cos(q(1))*cos(q(2))];

% Body frame to tail frame
R_bt = [cos(q(5)), 0, sin(q(5));
    sin(q(4))*sin(q(5)), cos(q(4)), -cos(q(5))*sin(q(4));
    -cos(q(4))*sin(q(5)), sin(q(4)), cos(q(4))*cos(q(5))];

% Inertial frame of reference to tail frame
R_Ot = R_Ob * R_bt;

%% Kinematics

% Body angular velocity
w_Ob__b = skew2angvel(R_Ob'*reshape(jtimes(reshape(R_Ob, 9, 1), q(1:3), q_dot(1:3)), 3, 3));

% Jacobian matrix of body angular velocity
J_Ob__b = jacobian(w_Ob__b, q_dot(1:3));

% Tail angular velocity
w_bt__b = skew2angvel(R_bt'*reshape(jacobian(reshape(R_bt, 9, 1), q(4:5))*q_dot(4:5), 3, 3));

% Jacobian matrix of tail angular velocity
J_bt__b = jacobian(w_bt__b, q_dot(4:5));

% Complete Jacobian
J__b = blkdiag(J_Ob__b, J_bt__b);

% Jacobian derivative
J_dot__b = reshape(jacobian(reshape(J__b, numel(J__b), 1), q)*q_dot, size(J__b, 1), size(J__b, 2));

%% Mass translation velocities

% COM to body and tail should be consistant with joint to body and tail
tmp = [diag(- ones(1, 3)), diag(ones(1, 3));
    diag(m_b .* ones(1, 3)), diag(m_t .* ones(1, 3))] \ ...
    [R_Ob * l_b - R_Ob * R_bt * l_t; zeros(3, 1)];

% Mass position w.r.t. COM
p_Ob__s = tmp(1:3, :);
p_Ot__s = tmp(4:6, :);

% Mass translation velocity w.r.t. COM
p_dot_Ob__s = jacobian(p_Ob__s, q) * q_dot;
p_dot_Ot__s = jacobian(p_Ot__s, q) * q_dot;

%% Lagrangian

% Body twist
V_Ob__b = [R_Ob'*p_dot_Ob__s; w_Ob__b];

% Tail twist, we need to add the body angular velocity here to get the full
% kinematics cause tail's rotation is on the body
V_Ot__b = [R_Ot'*p_dot_Ot__s; R_bt' * w_Ob__b + w_bt__b];

% Mass and inertia
M_b = diag([m_b * ones(3, 1); J_b]);
M_t = diag([m_t * ones(3, 1); J_t]);

% It's free falling or floating, so we don't need to consider gravity
T = V_Ob__b'*M_b*V_Ob__b/2+V_Ot__b' * M_t * V_Ot__b / 2;
V = 0;
L = T - V;

%% EOM

% Input in generalized coordinates, body should be zeros and tail's motor
% input should be same as roll and pitch
u = SX.sym('u', [5, 1]);

dLdq_dot = jacobian(L, q_dot);
dLdq = jacobian(L, q);

% Mass matrix
M = jacobian(dLdq_dot, q_dot);

% Other terms
h = jacobian(dLdq_dot, q)*q_dot - dLdq';

% Generalized coordinates acceleration
q_ddot = M\(u-h);

% Angular velocity change rate
w_dot = J__b*q_ddot+J_dot__b*q_dot;

% Function handle
dynamics_func = Function('dynamics_func',{q, q_dot, u},{w_dot},{'q','q_dot','u'},{'w_dot'});

%% Coordinates transformation

% Angular velocity of body w.r.t. inertial frame of reference in body frame
% and tail w.r.t. body frame in tail frame
w = SX.sym('w', [6, 1]);

% Map it back to Euler angle rates, singularity happens here
q_dot = pinv(J__b)*w;

w_dot = dynamics_func('q', q, 'q_dot', q_dot, 'u', u);
w_dot = w_dot.w_dot;

%% Continuous dynamics

% System states
x = [q; w];

% System dynamics
x_dot = [q_dot; w_dot];

dynamics_func = Function('dynamics_func',{x, u},{x_dot},{'x','u'},{'x_dot'});

% Dynamics jacobian
A = jacobian(x_dot, x);
B = jacobian(x_dot, u);

dynamics_jacobian_func = Function('dynamics_jacobian_func',{x, u},{A, B},{'x','u'},{'A', 'B'});

%% Discrete dynamics

% Discrete system states
xk = SX.sym('xk', [11, 1]);

% Discrete system inputs
uk = SX.sym('uk', [5, 1]);

% Discrete time step
dt = SX.sym('dt');

% RK4
k1 = dynamics_func('x', xk, 'u', uk);
k2 = dynamics_func('x', xk+dt*k1.x_dot/2, 'u', uk);
k3 = dynamics_func('x', xk+dt*k2.x_dot/2, 'u', uk);
k4 = dynamics_func('x', xk+dt*k3.x_dot, 'u', uk);

xk1 = xk + dt*(k1.x_dot + 2*k2.x_dot + 2*k3.x_dot + k4.x_dot)/6;

discrete_dynamics_func = Function('discrete_dynamics_func',{xk, uk, dt},{xk1},{'x','u','dt'},{'xk1'});

% Discrete dynamics jacobian
Ak = jacobian(xk1, xk);
Bk = jacobian(xk1, uk);

discrete_dynamics_jacobian_func = Function('discrete_dynamics_jacobian_func',{xk, uk, dt},{Ak, Bk},{'x','u','dt'},{'Ak', 'Bk'});

%% Output

model.dynamics = dynamics_func;
model.dynamics_jacobian = dynamics_jacobian_func;
model.discrete_dynamics = discrete_dynamics_func;
model.discrete_dynamics_jacobian = discrete_dynamics_jacobian_func;

end

