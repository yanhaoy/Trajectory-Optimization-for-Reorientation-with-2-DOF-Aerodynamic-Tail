clear; clc; close all

% state and input sizes
vals.n = 14; % full state space
vals.n_r = 12; % reduced quaternion state space
vals.m = 2; % inputs space

% set time span and knot points
vals.t0 = 0;
vals.tf = 5;
vals.N = 100;
vals.dt = (vals.tf-vals.t0) / vals.N; % nominal time step

vals.Nh = vals.n*vals.N + vals.m*(vals.N-1) + vals.N-1; % number of primals
[x_ind, u_ind, t_ind] = get_inds(vals); % get indices from packed states

% reference trajectory
vals.xref = zeros(vals.Nh,1); 

R0 = pi/4*rand(1,3) - pi/2; % random initial orientation
Rf = [0,0,0]; % final orientation
for i = 1:vals.N
    vals.xref(x_ind{i}(1:4)) = eul2quat((vals.N-i)/(vals.N-1)*R0+(i-1)/(vals.N-1)*Rf)';
    vals.xref(x_ind{i}(5:8)) = eul2quat(Rf)';
    
    if i < vals.N
        vals.xref(t_ind{i}) = vals.dt;
    end
end

% using inertial tail with quaternions
addpath('quaternion')
model_inert = init_dynamics();

% solve constrained optimization
options = optimoptions('fmincon','SpecifyConstraintGradient',true,...
    'SpecifyObjectiveGradient',true,'Display','iter','MaxFunctionEvaluations',250);
Z_inert = fmincon(@(x)cost_fun(x,vals),vals.xref,[],[],[],[],[],[],@(x)nonlcon(x,model_inert,vals),options);

% using heavy aerodynamic tail
addpath('quaternion_aerodynamics_ignore_translation')
model_heavy = init_dynamics_heavy();

% solve constrained optimization
options = optimoptions('fmincon','SpecifyConstraintGradient',true,...
    'SpecifyObjectiveGradient',true,'Display','iter','MaxFunctionEvaluations',250);
Z_heavy = fmincon(@(x)cost_fun(x,vals),vals.xref,[],[],[],[],[],[],@(x)nonlcon(x,model_heavy,vals),options);

% using light aerodynamic tail
addpath('quaternion_aerodynamics_ignore_translation')
model_light = init_dynamics_light();

% solve constrained optimization
options = optimoptions('fmincon','SpecifyConstraintGradient',true,...
    'SpecifyObjectiveGradient',true,'Display','iter','MaxFunctionEvaluations',250);
Z_light = fmincon(@(x)cost_fun(x,vals),vals.xref,[],[],[],[],[],[],@(x)nonlcon(x,model_light,vals),options);

%[cineq,ceq,Gineq,Geq]=nonlcon(Z,model,vals);
%viz_tail(Z,x_ind,u_ind,vals);
figure('Name','Intertial Tail States');
plot_states_comparison(Z_inert,Z_heavy,Z_light,x_ind,u_ind,t_ind,vals);

% cost function
function [f,g] = cost_fun(Z,vals)

[x_ind, u_ind, t_ind] = get_inds(vals);

f = 0; % cost function
g = zeros(vals.Nh,1); % gradient

% weight matrices
w = 1;
Q = eye(6);
R = 2*eye(vals.m);
S = 5;

for i = 1:vals.N
    q_b = Z(x_ind{i}(1:4)); % body quaternion
    q_t = Z(x_ind{i}(5:8)); % tail quaternion
    
    q0_b = vals.xref(x_ind{i}(1:4)); % reference body quaternion
    q0_t = vals.xref(x_ind{i}(5:8)); % reference tail quaternion
    
    q_b_cost = w*min(1-q_b'*q0_b,1+q_b'*q0_b); % simplified quaternion error 
    q_t_cost = w*min(1-q_t'*q0_t,1+q_t'*q0_t);
    
    % rotational velocity costs
    vel_cost = 0.5*(Z(x_ind{i}(9:end))-vals.xref(x_ind{i}(9:end)))'*Q*(Z(x_ind{i}(9:end))-vals.xref(x_ind{i}(9:end)));

    % state error gradients
    g(x_ind{i}(1:4)) = -w*sign(q_b'*q0_b)*q0_b;
    g(x_ind{i}(5:8)) = -w*sign(q_t'*q0_t)*q0_t;
    g(x_ind{i}(9:end)) = Q*(Z(x_ind{i}(9:end))-vals.xref(x_ind{i}(9:end)));
    
    if i < vals.N
        % input cost and gradient
        u_cost = 0.5*(Z(u_ind{i}) - vals.xref(u_ind{i}))'*R*(Z(u_ind{i}) - vals.xref(u_ind{i}));
        g(u_ind{i}) = R*(Z(u_ind{i}) - vals.xref(u_ind{i}));
        
        % minimum time cost and gradient
        t_cost = S*Z(t_ind{i});
        g(t_ind{i}) = S;
    else
        u_cost = 0;
        t_cost = 0;
    end
    
    f = f + q_b_cost + q_t_cost + vel_cost + u_cost + t_cost; % full cost
end

end

% constraints
function [cineq,ceq,Gineq,Geq]=nonlcon(Z,model,vals)

n = vals.n; m = vals.m; N = vals.N; n_r = vals.n_r;
[x_ind, u_ind, t_ind] = get_inds(vals);

ceq = zeros(n_r*(N+1)+2*N,1);
Geq = zeros(vals.Nh,n_r*(N+1)+2*N);
cineq = zeros(N-1,1);
Gineq = zeros(vals.Nh,N-1);

for i = 1:N-1
    % normalize quaternions
    x = Z(x_ind{i}); x(1:4) = x(1:4)/norm(x(1:4)); 
    x(5:8) = x(5:8)/norm(x(5:8));
    
    u = Z(u_ind{i});
    
    x_dyn = discrete_dynamics(model,x,u,Z(t_ind{i}));
    
    % dynamics contraints
    ceq(n_r*(i-1)+1:n_r*(i-1)+n_r) = state_error(x_dyn, Z(x_ind{i+1}));
    
    % unit quaternion constraints
    ceq(n_r*(N+1)+i) = norm(Z(x_ind{i}(1:4))) - 1;
    ceq(n_r*(N+1)+N+i) = norm(Z(x_ind{i}(5:8))) - 1;
    
    % positive time step constraints
    cineq(i) = -Z(t_ind{i});
    Gineq(t_ind{i},i) = -1;
    
    % constraint Jacobians
    [Ak, Bk] = discrete_dynamics_jacobian(model, x, u, Z(t_ind{i}));
    E = state_error_jacobian(x,vals);
    
    Ak_r = E'*Ak;
    Bk_r = E'*Bk;
    
    Geq(x_ind{i},n_r*(i-1)+1:n_r*(i-1)+n_r) = Ak_r';
    Geq(u_ind{i},n_r*(i-1)+1:n_r*(i-1)+n_r) = Bk_r';
    Geq(x_ind{i+1},n_r*(i-1)+1:n_r*(i-1)+n_r) = -E;
    
    Geq(x_ind{i}(1:4),n_r*(N+1)+i) = Z(x_ind{i}(1:4)) ./ norm(Z(x_ind{i}(1:4)));
    Geq(x_ind{i}(5:8),n_r*(N+1)+N+i) = Z(x_ind{i}(5:8)) ./ norm(Z(x_ind{i}(5:8)));

end

% initial state constraint
ceq(n_r*(N-1)+1:n_r*(N-1)+n_r) = state_error(Z(x_ind{1}), vals.xref(1:n));
E = state_error_jacobian(Z(x_ind{1}),vals);
Geq(x_ind{1},n_r*(N-1)+1:n_r*(N-1)+n_r) = E;

% final state constraint
ceq(n_r*N+1:n_r*N+n_r) = state_error(Z(x_ind{end}), vals.xref(x_ind{end}));
E = state_error_jacobian(Z(x_ind{end}),vals);
Geq(x_ind{end},n_r*N+1:n_r*N+n_r) = E;

% final unit quaternion constraint
ceq(n_r*(N+1)+N) = norm(Z(x_ind{end}(1:4))) - 1;
ceq(n_r*(N+1)+2*N) = norm(Z(x_ind{end}(5:8))) - 1;

Geq(x_ind{end}(1:4),n_r*(N+1)+N) = Z(x_ind{end}(1:4)) ./ norm(Z(x_ind{end}(1:4)));
Geq(x_ind{end}(5:8),n_r*(N+1)+2*N) = Z(x_ind{end}(5:8)) ./ norm(Z(x_ind{end}(5:8)));
end

% get indices for states, inputs, and time steps
function [x_ind, u_ind, t_ind] = get_inds(vals)

n = vals.n; m = vals.m; N = vals.N; n_r = vals.n_r;

x_ind = cell(N,1);
u_ind = cell(N-1,1);
t_ind = cell(N-1,1);

for i = 1:N-1
    x_ind{i} = 1+(i-1)*(n+m+1):n+(i-1)*(n+m+1);
    u_ind{i} = n+1+(i-1)*(n+m+1):n+m+(i-1)*(n+m+1);
    t_ind{i} = n+m+1+(i-1)*(n+m+1);
end
x_ind{N} = 1+(N-1)*(n+m+1):n+(N-1)*(n+m+1);

end

% discrete dynamics
function xk1 = discrete_dynamics(model, x, u, dt)

tmp = model.discrete_dynamics('x', x, 'u', u, 'dt', dt);
xk1 = full(tmp.xk1);

end

% discrete dynamics Jacobian
function [Ak, Bk] = discrete_dynamics_jacobian(model, x, u, dt)

tmp = model.discrete_dynamics_jacobian('x', x, 'u', u, 'dt', dt);
Ak = full(tmp.Ak);
Bk = full(tmp.Bk);

end

% operator for quaternion error
function L = lmult(q)

v_hat = [0, -q(4), q(3); q(4), 0, -q(2); -q(3), q(2), 0];
L = [q(1), -q(2:4)'; q(2:4), q(1)*eye(3)+v_hat];

end

% calculate state error with 3-parameter expression for quaternions
function dx = state_error(x,x0)

dx = zeros(length(x)-2,1);
q_diff_b = lmult(x0(1:4))'*x(1:4);
q_diff_t = lmult(x0(5:8))'*x(5:8);

dx(1:3) = 1/q_diff_b(1)*q_diff_b(2:4);
dx(4:6) = 1/q_diff_t(1)*q_diff_t(2:4);
dx(7:end) = x(9:end) - x0(9:end);

end

% state error Jacobian to reduced quaternion coordinates
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