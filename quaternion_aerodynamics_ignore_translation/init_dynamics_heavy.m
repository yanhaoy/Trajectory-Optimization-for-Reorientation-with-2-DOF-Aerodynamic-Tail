function [model] = init_dynamics()

% addpath('..\utils\')
% addpath('..\utils\casadi\')
addpath('utils')
addpath('utils\casadi')
import casadi.*

% Body mass
m_b = 1;

% Tail mass, since we are using aerodyanmics, we could have a lighter tail
% m_t = 0.25;
m_t = 0.1;

% Body inertia
J_b = [1; 1; 1];

% Tail inertia, since we are using aerodyanmics, we could have a lighter
% tail
% J_t = [0.25; 0.25; 0.25];
J_t = [0.1;0.1;0.1];

% Translation from body to joint in body frame
l_b = [1; 0; 0];

% Tail length, since we are using aerodyanmics, we could have a longer tail
tail_l = 1;

% Just a trick to get the relationship between aerodynamics magnitude and
% length for the integration
length = SX.sym('length');
l_t = [- length; 0; 0];

% Air density
rho = 1.29;

% Drag coefficient that captures the surface interaction between the fluid 
% and the body
C_D = 2;

% Tail shaft width
shaft_w = 1;

% Length from joint to tail shaft
shaft_l = 0;

% Generalized coordinates: quaternions from inertial frame of reference to
% body; quaternions from body to tail
q = SX.sym('q', [8, 1]);

% Quaternions change rates
q_dot = SX.sym('q_dot', [8, 1]);

%% Rotation matrices

% Inertial frame of reference to body frame
R_Ob = quat2rotm(q(1:4));

% Body frame to tail frame
R_bt = quat2rotm(q(5:8));

% Inertial frame of reference to tail frame
R_Ot = R_Ob * R_bt;

%% Kinematics

% Body angular velocity
w_Ob__b = 2*quat2L(q(1:4))*q_dot(1:4);

% Tail angular velocity
w_bt__b = 2*quat2L(q(5:8))*q_dot(5:8);

%% Mass translation velocities

% COM to body and tail should be consistant with joint to body and tail
tmp = [diag(- ones(1, 3)), diag(ones(1, 3));
    diag(m_b .* ones(1, 3)), diag(m_t .* ones(1, 3))] \ ...
    [R_Ob * l_b - R_Ob * R_bt * l_t; zeros(3, 1)];

% Mass position w.r.t. COM
p_Ob__s = tmp(1:3, :);
p_Ot__s = tmp(4:6, :);

% Mass translation velocity w.r.t. COM
p_dot_Ob__s = jtimes(p_Ob__s, q, q_dot);
p_dot_Ot__s = jtimes(p_Ot__s, q, q_dot);

%% Aerodynamics

% Here we ignore the translational effect like the Aerodynamics tail paper, 
% because the direction of the external force here keeps changing, and the 
% integrated speed should tend to zero.

% We know the translation velocity is a linear function of the length
k_a = jacobian(p_dot_Ot__s, length);
k_b = p_dot_Ot__s - k_a*length;

% This is the integration result along the shaft
magnitude = rho*C_D*shaft_w/2* ...
    (k_a'*k_a*(tail_l^4 - shaft_l^4)/4 + ...
    2*k_a'*k_b*(tail_l^3 - shaft_l^3)/3 + ...
    k_b'*k_b*(tail_l^2 - shaft_l^2)/2);

% Maybe there's a better way to compute both magnitude and direction
% simultaneously, but currently I have it seperate. Here I have to add a
% small constant to avoid numerical issues
direction = cross(l_t/sqrt(l_t'*l_t + 1e-12), ...
    R_Ot'*p_dot_Ot__s/sqrt(p_dot_Ot__s'*p_dot_Ot__s + 1e-12));

% Resulted aerodynamic torque in tail quaternion coordinates
u_ad_t = 2*quat2L(q(5:8))'*direction*magnitude;

% Resulted aerodynamic torque in body quaternion coordinates
u_ad_b = 2*quat2L(q(1:4))'*R_bt*direction*magnitude;

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

% Input in generalized coordinates
u = SX.sym('u', [8, 1]);

dLdq_dot = jacobian(L, q_dot);
dLdq = jacobian(L, q);

% Mass matrix
M = jacobian(dLdq_dot, q_dot);

% Other terms
h = jacobian(dLdq_dot, q)*q_dot - dLdq';

% Constraints Jacobian
J_c = blkdiag(2*q(1:4)', 2*q(5:8)');
h_c = [2*q_dot(1:4)'*q_dot(1:4); 2*q_dot(5:8)'*q_dot(5:8)];

% Generalized coordinates acceleration
q_ddot = [M, J_c'; J_c, zeros(2, 2)]\([u-h; zeros(2, 1)]-[zeros(8, 1); h_c]+[u_ad_b; u_ad_t; zeros(2, 1)]);
q_ddot = q_ddot(1:8);

% Angular velocity change rate
w_dot = blkdiag(2*quat2L(q(1:4)), 2*quat2L(q(5:8)))*q_ddot;

% Function handle
dynamics_func = Function('dynamics_func',{q, q_dot, u, length},{w_dot},{'q','q_dot','u','length'},{'w_dot'});

%% Coordinates transformation

% Angular velocity of body w.r.t. inertial frame of reference in body frame
% and tail w.r.t. body frame in tail frame
w = SX.sym('w', [6, 1]);

% Input of tail's roll and pitch motor
u = SX.sym('u', [2, 1]);

% Pitch motor position
pitch = atan2(R_bt(1, 3), R_bt(1, 1));

% Jacobian mapping u to torque in body frame
J_u = [cos(pitch), 0;
    0, 1;
    sin(pitch), 0];

% Jacobian mapping u to quaternion coordinates
J_u = 2*quat2L(q(5:8))'*J_u;

% Map it back to quaternions rates, no singularity here
q_dot = blkdiag(quat2L(q(1:4))'/2, quat2L(q(5:8))'/2)*w;

w_dot = dynamics_func('q', q, 'q_dot', q_dot, 'u', [zeros(4, 1); J_u*u], 'length', tail_l);
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
xk = SX.sym('xk', [14, 1]);

% Discrete system inputs
uk = SX.sym('uk', [2, 1]);

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
Tk = jacobian(xk1, dt);

discrete_dynamics_jacobian_func = Function('discrete_dynamics_jacobian_func',{xk, uk, dt},{Ak, Bk, Tk},{'x','u','dt'},{'Ak', 'Bk', 'Tk'});

%% Output

model.dynamics = dynamics_func;
model.dynamics_jacobian = dynamics_jacobian_func;
model.discrete_dynamics = discrete_dynamics_func;
model.discrete_dynamics_jacobian = discrete_dynamics_jacobian_func;

end

