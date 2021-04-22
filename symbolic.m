clear; clc; close all
addpath('utils')

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
% body, roll, pitch, and yaw; Euler angles from body to tail, roll, pitch,
% and yaw
syms q [6, 1] real

% Angular velocities: body angular velocities w.r.t. inertial frame of
% reference, represented in body frame; tail angular velocities w.r.t. body
% frame, represented in tail frame
syms w [6, 1] real

%% Rotation matrices

R_roll_Ob = expm(angvel2skew([1; 0; 0]) * q(1));
R_pitch_Ob = expm(angvel2skew([0; 1; 0]) * q(2));
R_yaw_Ob = expm(angvel2skew([0; 0; 1]) * q(3));

% Inertial frame of reference to body frame
R_Ob = R_roll_Ob * R_pitch_Ob * R_yaw_Ob;

R_roll_bt = expm(angvel2skew([1; 0; 0]) * q(4));
R_pitch_bt = expm(angvel2skew([0; 1; 0]) * q(5));
R_yaw_bt = expm(angvel2skew([0; 0; 1]) * q(6));

% body frame to tail frame
R_bt = R_roll_bt * R_pitch_bt * R_yaw_bt;

% Inertial frame of reference to tail frame
R_Ot = R_Ob * R_bt;

%% Jacobian matrices

% Jacobian matrix of body angular velocity
J_Ob__b = jacobian(reshape(R_Ob, [], 1), q(1:3));
J_Ob__b = [skew2angvel(R_Ob'*reshape(J_Ob__b(:, 1), 3, 3)), ...
    skew2angvel(R_Ob'*reshape(J_Ob__b(:, 2), 3, 3)), ...
    skew2angvel(R_Ob'*reshape(J_Ob__b(:, 3), 3, 3))];

% Jacobian matrix of tail angular velocity
J_bt__b = jacobian(reshape(R_bt, [], 1), q(4:6));
J_bt__b = [skew2angvel(R_bt'*reshape(J_bt__b(:, 1), 3, 3)), ...
    skew2angvel(R_bt'*reshape(J_bt__b(:, 2), 3, 3)), ...
    skew2angvel(R_bt'*reshape(J_bt__b(:, 3), 3, 3))];

% Complete Jacobian
J__b = simplify(blkdiag(J_Ob__b, J_bt__b));

% Equivalent Euler angle change rate
q_dot = J__b \ w;

%% Mass translation velocities

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
V_Ob__b = simplify([R_Ob'*p_dot_Ob__s; w(1:3)]);

% Tail twist, we need to add the body angular velocity here to get the full
% kinematics cause tail's rotation is on the body
V_Ot__b = simplify([R_Ot'*p_dot_Ot__s; R_bt' * w(1:3) + w(4:6)]);

% Mass and inertia
M_b = diag([m_b * ones(3, 1); J_b]);
M_t = diag([m_t * ones(3, 1); J_t]);

% It's free falling or floating, so we don't need to consider gravity
T = V_Ob__b'*M_b*V_Ob__b/2+V_Ot__b' * M_t * V_Ot__b / 2;
V = 0;
L = T - V;

%% EOM

% Since we have Euler angle as the generalized coordinates but use angular
% velocity instead of the Euler angle rate, we have to use the chain rule
% to compute the COM.

% d(L) / d(q_dot) = d(L) / d(w) * d(w) / d(q_dot)
dLdq_dot = simplify(jacobian(L, w) * J__b);

% d(L(w(q))) / d(q) = d(L(w)) / d(q) + d(L(w)) / d(w) * d(w) / d(q)
dwdq = jacobian(reshape(J__b, [], 1), q);
dwdq = [reshape(dwdq(:, 1), 6, 6) * q_dot, ...
    reshape(dwdq(:, 2), 6, 6) * q_dot, ...
    reshape(dwdq(:, 3), 6, 6) * q_dot, ...
    reshape(dwdq(:, 4), 6, 6) * q_dot, ...
    reshape(dwdq(:, 5), 6, 6) * q_dot, ...
    reshape(dwdq(:, 6), 6, 6) * q_dot];
dLdq = simplify(jacobian(L, q) + jacobian(L, w) * dwdq);

% d(d(L) / d(q_dot)) / d(t) - d(L) / d(q)
% = d(d(L) /d(q_dot)) / d(w) * w_dot + d(q_dot)) / d(q) * q_dot - d(L) / d(q)
% = M * w_dot + h
M = simplify(jacobian(dLdq_dot, w));
h = simplify(jacobian(dLdq_dot, q) * q_dot - dLdq');

%%
matlabFunction(M, h, q_dot, 'File', 'compute_dynamics', 'Vars', {q, w});

%%
dM = simplify(jacobian(reshape(M, [], 1), [q; w]));
dh = simplify(jacobian(reshape(h, [], 1), [q; w]));

%%
matlabFunction(dM, dh, 'File', 'compute_dynamics_jacobian', 'Vars', {q, w});