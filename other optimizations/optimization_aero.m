clear; clc; close all

% reorienting body with inertial tail
%addpath('utils')
%addpath('utils\casadi')
%import casadi.*

% state and input sizes
vals.n = 14;
vals.n_r = 12;
vals.m = 2;

% set time span and knot points
vals.t0 = 0;
vals.tf = 5;
vals.N = 100;
vals.dt = (vals.tf-vals.t0) / vals.N;
vals.Nh = vals.n*vals.N + vals.m*(vals.N-1);
[x_ind, u_ind] = get_inds(vals);

% initial and final conditions
vals.x0 = zeros(vals.Nh,1); 
vals.xf = zeros(vals.Nh,1); 

R0 = pi/4*rand(1,3) - pi/2;
Rf = [0,0,0];
for i = 1:vals.N
    vals.x0(1+(i-1)*(vals.n+vals.m):4+(i-1)*(vals.n+vals.m)) = eul2quat(R0)';
    vals.x0(5+(i-1)*(vals.n+vals.m):8+(i-1)*(vals.n+vals.m)) = eul2quat(Rf)';
    vals.xf(1+(i-1)*(vals.n+vals.m):4+(i-1)*(vals.n+vals.m)) = eul2quat((vals.N-i)/(vals.N-1)*R0+(i-1)/(vals.N-1)*Rf)';
    vals.xf(5+(i-1)*(vals.n+vals.m):8+(i-1)*(vals.n+vals.m)) = eul2quat(Rf)';
end

% using quaternions
addpath('quaternion_aerodynamics_ignore_translation')
model = init_dynamics();

% solve constrained optimization
options = optimoptions('fmincon','SpecifyConstraintGradient',true,...
    'SpecifyObjectiveGradient',true,'Display','iter','MaxFunctionEvaluations',250);
Z = fmincon(@(x)cost_fun(x,vals),vals.xf,[],[],[],[],[],[],@(x)nonlcon(x,model,vals),options);

[cineq,ceq,Gineq,Geq]=nonlcon(Z,model,vals);
%viz_tail(Z,x_ind,u_ind,vals);
plot_states(Z,x_ind,u_ind,vals);

% cost function
function [f,g] = cost_fun(Z,vals)

[x_ind, u_ind] = get_inds(vals);

f = 0;
g = zeros(vals.Nh,1);

% weight matrices
w = 1;
P = eye(6);
R = 2*eye(vals.m);

for i = 1:vals.N
    q_b = Z(x_ind{i}(1:4));
    q_t = Z(x_ind{i}(5:8));
    
    q0_b = vals.xf(x_ind{i}(1:4));
    q0_t = vals.xf(x_ind{i}(5:8));
    
    q_b_cost = w*min(1-q_b'*q0_b,1+q_b'*q0_b);
    q_t_cost = w*min(1-q_t'*q0_t,1+q_t'*q0_t);
    vel_cost = 0.5*(Z(x_ind{i}(9:end))-vals.xf(x_ind{i}(9:end)))'*P*(Z(x_ind{i}(9:end))-vals.xf(x_ind{i}(9:end)));

    g(x_ind{i}(1:4)) = -w*sign(q_b'*q0_b)*q0_b;
    g(x_ind{i}(5:8)) = -w*sign(q_t'*q0_t)*q0_t;
    g(x_ind{i}(9:end)) = P*(Z(x_ind{i}(9:end))-vals.xf(x_ind{i}(9:end)));
    
    if i < vals.N
        u_cost = 0.5*(Z(u_ind{i}) - vals.xf(u_ind{i}))'*R*(Z(u_ind{i}) - vals.xf(u_ind{i}));
        g(u_ind{i}) = R*(Z(u_ind{i}) - vals.xf(u_ind{i}));
    else
        u_cost = 0;
    end
    
    f = f + q_b_cost + q_t_cost + vel_cost + u_cost;
end

end

% constaints
function [cineq,ceq,Gineq,Geq]=nonlcon(Z,model,vals)

n = vals.n; m = vals.m; N = vals.N; n_r = vals.n_r;
[x_ind, u_ind, x_ind_r, u_ind_r] = get_inds(vals);

ceq = zeros(n_r*(N+1)+2*N,1);
Geq = zeros(vals.Nh,n_r*(N+1)+2*N);
cineq = [];
Gineq = [];

for i = 1:N-1
    x = Z(x_ind{i}); x(1:4) = x(1:4)/norm(x(1:4)); x(5:8) = x(5:8)/norm(x(5:8));
    u = Z(u_ind{i});
    
    x_dyn = discrete_dynamics(model,x,u,vals.dt);
    
    % dynamics contraints
    ceq(n_r*(i-1)+1:n_r*(i-1)+n_r) = state_error(x_dyn, Z(x_ind{i+1}));
    
    % unit quaternion constraints
    ceq(n_r*(N+1)+i) = norm(Z(x_ind{i}(1:4))) - 1;
    ceq(n_r*(N+1)+N+i) = norm(Z(x_ind{i}(5:8))) - 1;
    
    [Ak, Bk] = discrete_dynamics_jacobian(model, x, u, vals.dt);
    E = state_error_jacobian(x,vals);
    
    Ak_r = E'*Ak;
    Bk_r = E'*Bk;
    
    Geq(x_ind{i},n_r*(i-1)+1:n_r*(i-1)+n_r) = Ak_r';
    Geq(u_ind{i},n_r*(i-1)+1:n_r*(i-1)+n_r) = Bk_r';
    Geq(x_ind{i+1},n_r*(i-1)+1:n_r*(i-1)+n_r) = -E;
    
    Geq(x_ind{i}(1:4),n_r*(N+1)+i) = Z(x_ind{i}(1:4)) ./ norm(Z(x_ind{i}(1:4)));
    Geq(x_ind{i}(5:8),n_r*(N+1)+N+i) = Z(x_ind{i}(5:8)) ./ norm(Z(x_ind{i}(5:8)));

end

ceq(n_r*(N-1)+1:n_r*(N-1)+n_r) = state_error(Z(x_ind{1}), vals.xf(1:n));
E = state_error_jacobian(Z(x_ind{1}),vals);
Geq(x_ind{1},n_r*(N-1)+1:n_r*(N-1)+n_r) = E;

ceq(n_r*N+1:n_r*N+n_r) = state_error(Z(x_ind{end}), vals.xf(end-n+1:end));
E = state_error_jacobian(Z(x_ind{end}),vals);
Geq(x_ind{end},n_r*N+1:n_r*N+n_r) = E;

ceq(n_r*(N+1)+N) = norm(Z(x_ind{end}(1:4))) - 1;
ceq(n_r*(N+1)+2*N) = norm(Z(x_ind{end}(5:8))) - 1;

Geq(x_ind{end}(1:4),n_r*(N+1)+N) = Z(x_ind{end}(1:4)) ./ norm(Z(x_ind{end}(1:4)));
Geq(x_ind{end}(5:8),n_r*(N+1)+2*N) = Z(x_ind{end}(5:8)) ./ norm(Z(x_ind{end}(5:8)));
end

function [x_ind, u_ind, x_ind_r, u_ind_r] = get_inds(vals)

n = vals.n; m = vals.m; N = vals.N; n_r = vals.n_r;

x_ind = cell(N,1);
u_ind = cell(N,1);
x_ind_r = cell(N,1);
u_ind_r = cell(N,1);

for i = 1:N-1
    x_ind{i} = 1+(i-1)*(n+m):n+(i-1)*(n+m);
    u_ind{i} = n+1+(i-1)*(n+m):n+m+(i-1)*(n+m);
    x_ind_r{i} = 1+(i-1)*(n_r+m):n_r+(i-1)*(n_r+m);
    u_ind_r{i} = n_r+1+(i-1)*(n_r+m):n_r+m+(i-1)*(n_r+m);
end
x_ind{N} = 1+(N-1)*(n+m):n+(N-1)*(n+m);
x_ind_r{N} = 1+(N-1)*(n_r+m):n_r+(N-1)*(n_r+m);
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

function L = lmult(q)

v_hat = [0, -q(4), q(3); q(4), 0, -q(2); -q(3), q(2), 0];
L = [q(1), -q(2:4)'; q(2:4), q(1)*eye(3)+v_hat];

end

function dx = state_error(x,x0)

dx = zeros(length(x)-2,1);
q_diff_b = lmult(x0(1:4))'*x(1:4);
q_diff_t = lmult(x0(5:8))'*x(5:8);

dx(1:3) = 1/q_diff_b(1)*q_diff_b(2:4);
dx(4:6) = 1/q_diff_t(1)*q_diff_t(2:4);
dx(7:end) = x(9:end) - x0(9:end);

end

function J = state_error_jacobian(x,vals)

n = vals.n; m = vals.m; N = vals.N; n_r = vals.n_r;

J = zeros(n,n_r);

q_b = x(1:4);
v_hat_b = [0, -q_b(4), q_b(3); q_b(4), 0, -q_b(2); -q_b(3), q_b(2), 0];
G_b = [-q_b(2:4)'; q_b(1)*eye(3) + v_hat_b];

q_t = x(5:8);
v_hat_t = [0, -q_t(4), q_t(3); q_t(4), 0, -q_t(2); -q_t(3), q_t(2), 0];
G_t = [-q_t(2:4)'; q_t(1)*eye(3) + v_hat_t];

J(1:4,1:3) = G_b;
J(5:8,4:6) = G_t;
J(9:n,7:n_r) = eye(6);

end