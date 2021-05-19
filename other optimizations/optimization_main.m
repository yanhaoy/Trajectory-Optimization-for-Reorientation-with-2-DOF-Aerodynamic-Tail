clear; clc; close all

% reference state and input sizes
vals.n{1} = 14; % full state space
vals.n{2} = 20; % full state space with body COM
vals.n_r = 12; % reduced quaternion state space
vals.m = 2; % input space

vals.mode = 1;

% set time span and knot points
vals.t0 = 0;
vals.tf = 10;
vals.N = 100;
vals.dt = (vals.tf-vals.t0) / vals.N; % nominal time step

% number of primals
vals.Nh{1} = vals.n{1}*vals.N + vals.m*(vals.N-1) + vals.N-1;
vals.Nh{2} = vals.n{2}*vals.N + vals.m*(vals.N-1) + vals.N-1;

% get indices from packed states
% first set mainly for plotting; second set used in optimization
[x_ind, u_ind, t_ind, x2_ind, u2_ind, t2_ind] = get_inds(vals);

% reference trajectory
vals.xref{1} = zeros(vals.Nh{1},1); 

R0 = [pi/6 pi/6 pi/6]; %pi/2*rand(1,3) - pi/4; % random initial orientation
Rf = [0,0,0]; % final orientation
for i = 1:vals.N
    vals.xref{1}(x_ind{i}(1:4)) = eul2quat((vals.N-i)/(vals.N-1)*R0+(i-1)/(vals.N-1)*Rf,'XYZ')';
    vals.xref{1}(x_ind{i}(5:8)) = eul2quat(Rf,'XYZ')';
    
    if i < vals.N
        vals.xref{1}(t_ind{i}) = vals.dt;
    end
end

% using inertial tail with quaternions
addpath('quaternion')
model_inert = init_dynamics();
vals.mode = 1;

% solve constrained optimization
options = optimoptions('fmincon','SpecifyConstraintGradient',true,...
    'SpecifyObjectiveGradient',true,'Display','iter','MaxFunctionEvaluations',250);
Z_inert = fmincon(@(x)cost_fun(x,vals),vals.xref{1},[],[],[],[],[],[],@(x)nonlcon(x,model_inert,vals),options);

% using heavy aerodynamic tail
addpath('quaternion_aerodynamics_free_falling')
model_heavy = init_dynamics_heavy();
vals.mode = 2;

% convert reference trajectory to larger state space
[x_ind, u_ind, t_ind, x2_ind, u2_ind, t2_ind] = get_inds(vals);
vals.xref{2} = zeros(vals.Nh{2},1); 

for i = 1:vals.N
    vals.xref{2}(x2_ind{i}(1:8)) = vals.xref{1}(x_ind{i}(1:8));
    vals.xref{2}(x2_ind{i}(12:17)) = vals.xref{1}(x_ind{i}(9:14));
    
    if i < vals.N
        vals.xref{2}(t2_ind{i}) = vals.xref{1}(t_ind{i});
    end
end

% solve constrained optimization
options = optimoptions('fmincon','SpecifyConstraintGradient',true,...
    'SpecifyObjectiveGradient',true,'Display','iter','MaxFunctionEvaluations',250);
Z_heavy = fmincon(@(x)cost_fun(x,vals),vals.xref{2},[],[],[],[],[],[],@(x)nonlcon(x,model_heavy,vals),options);

% using light aerodynamic tail
addpath('quaternion_aerodynamics_free_falling')
model_light = init_dynamics_light();
vals.mode = 2;

% solve constrained optimization
options = optimoptions('fmincon','SpecifyConstraintGradient',true,...
    'SpecifyObjectiveGradient',true,'Display','iter','MaxFunctionEvaluations',250);
Z_light = fmincon(@(x)cost_fun(x,vals),vals.xref{2},[],[],[],[],[],[],@(x)nonlcon(x,model_light,vals),options);

%[cineq,ceq,Gineq,Geq]=nonlcon(Z,model,vals);
%viz_tail(Z,x_ind,u_ind,vals);
figure('Name','Intertial Tail States');
plot_states_comparison(Z_inert,Z_heavy,Z_light,x_ind,u_ind,t_ind,x2_ind,u2_ind,t2_ind,vals);

% cost function
function [f,g] = cost_fun(Z,vals)

[~,~,~, x2_ind, u2_ind, t2_ind] = get_inds(vals);

f = 0; % cost function

g = zeros(vals.Nh{vals.mode},1); % gradient

% weight matrices
w = 2;
Q = 2*eye(6);
R = 3*eye(vals.m);
S = 10;

for i = 1:vals.N
    
    if i < vals.N
        
        h = S*Z(t2_ind{i});
        
        % input cost and gradient
        u_cost = h*0.5*(Z(u2_ind{i}) - vals.xref{vals.mode}(u2_ind{i}))'*R*(Z(u2_ind{i}) - vals.xref{vals.mode}(u2_ind{i}));
        g(u2_ind{i}) = h*R*(Z(u2_ind{i}) - vals.xref{vals.mode}(u2_ind{i}));
  
    else
        u_cost = 0;
        h = 1;
    end
    
    q_b = Z(x2_ind{i}(1:4)); % body quaternion
    q_t = Z(x2_ind{i}(5:8)); % tail quaternion
    
    q0_b = vals.xref{vals.mode}(x2_ind{i}(1:4)); % reference body quaternion
    q0_t = vals.xref{vals.mode}(x2_ind{i}(5:8)); % reference tail quaternion
    
    q_b_cost = h*w*min(1-q_b'*q0_b,1+q_b'*q0_b); % simplified quaternion error 
    q_t_cost = h*w*min(1-q_t'*q0_t,1+q_t'*q0_t);
    
    if vals.mode == 1
        rot_inds = 9:14;
    elseif vals.mode == 2
        rot_inds = 12:17;
    end
    
    % rotational velocity costs
    vel_cost = h*0.5*(Z(x2_ind{i}(rot_inds))-vals.xref{vals.mode}(x2_ind{i}(rot_inds)))'*Q*(Z(x2_ind{i}(rot_inds))-vals.xref{vals.mode}(x2_ind{i}(rot_inds)));

    % state error gradients
    g(x2_ind{i}(1:4)) = -h*w*sign(q_b'*q0_b)*q0_b;
    g(x2_ind{i}(5:8)) = -h*w*sign(q_t'*q0_t)*q0_t;
    g(x2_ind{i}(rot_inds)) = h*Q*(Z(x2_ind{i}(rot_inds))-vals.xref{vals.mode}(x2_ind{i}(rot_inds)));
    
    if i < vals.N
        g(t2_ind{i}) = (q_b_cost + q_t_cost + vel_cost + u_cost) / Z(t2_ind{i});
    end
    
    f = f + q_b_cost + q_t_cost + vel_cost + u_cost; % full cost
end

end

% constraints
function [cineq,ceq,Gineq,Geq]=nonlcon(Z,model,vals)

n = vals.n; m = vals.m; N = vals.N; n_r = vals.n_r;
[~,~,~, x2_ind, u2_ind, t2_ind] = get_inds(vals);

ceq = zeros(n_r*(N+1)+2*N,1);
Geq = zeros(vals.Nh{vals.mode},n_r*(N+1)+2*N);
cineq = zeros(N-1,1);
Gineq = zeros(vals.Nh{vals.mode},N-1);

for i = 1:N-1
    % normalize quaternions
    x = Z(x2_ind{i}); x(1:4) = x(1:4)/norm(x(1:4)); 
    x(5:8) = x(5:8)/norm(x(5:8));
    
    u = Z(u2_ind{i});
    
    x_dyn = discrete_dynamics(model,x,u,Z(t2_ind{i}));
    
    % dynamics contraints
    ceq(n_r*(i-1)+1:n_r*(i-1)+n_r) = state_error(x_dyn, Z(x2_ind{i+1}),vals);
    
    % unit quaternion constraints
    ceq(n_r*(N+1)+i) = norm(Z(x2_ind{i}(1:4))) - 1;
    ceq(n_r*(N+1)+N+i) = norm(Z(x2_ind{i}(5:8))) - 1;
    
    % positive time step constraints
    cineq(i) = -Z(t2_ind{i});
    Gineq(t2_ind{i},i) = -1;
    
    % constraint Jacobians
    [Ak, Bk, Tk] = discrete_dynamics_jacobian(model, x, u, Z(t2_ind{i}));
    
    E = state_error_jacobian(x,vals);
    
    Ak_r = E'*Ak;
    Bk_r = E'*Bk;
    Tk_r = E'*Tk;
    
    Geq(x2_ind{i},n_r*(i-1)+1:n_r*(i-1)+n_r) = Ak_r';
    Geq(u2_ind{i},n_r*(i-1)+1:n_r*(i-1)+n_r) = Bk_r';
    Geq(t2_ind{i},n_r*(i-1)+1:n_r*(i-1)+n_r) = Tk_r';
    Geq(x2_ind{i+1},n_r*(i-1)+1:n_r*(i-1)+n_r) = -E;
    
    Geq(x2_ind{i}(1:4),n_r*(N+1)+i) = Z(x2_ind{i}(1:4)) ./ norm(Z(x2_ind{i}(1:4)));
    Geq(x2_ind{i}(5:8),n_r*(N+1)+N+i) = Z(x2_ind{i}(5:8)) ./ norm(Z(x2_ind{i}(5:8)));

end

% initial state constraint
ceq(n_r*(N-1)+1:n_r*(N-1)+n_r) = state_error(Z(x2_ind{1}), vals.xref{vals.mode}(x2_ind{1}),vals);
E = state_error_jacobian(Z(x2_ind{1}),vals);
Geq(x2_ind{1},n_r*(N-1)+1:n_r*(N-1)+n_r) = E;

% final state constraint
ceq(n_r*N+1:n_r*N+n_r) = state_error(Z(x2_ind{end}), vals.xref{vals.mode}(x2_ind{end}),vals);
E = state_error_jacobian(Z(x2_ind{end}),vals);
Geq(x2_ind{end},n_r*N+1:n_r*N+n_r) = E;

% final unit quaternion constraint
ceq(n_r*(N+1)+N) = norm(Z(x2_ind{end}(1:4))) - 1;
ceq(n_r*(N+1)+2*N) = norm(Z(x2_ind{end}(5:8))) - 1;

Geq(x2_ind{end}(1:4),n_r*(N+1)+N) = Z(x2_ind{end}(1:4)) ./ norm(Z(x2_ind{end}(1:4)));
Geq(x2_ind{end}(5:8),n_r*(N+1)+2*N) = Z(x2_ind{end}(5:8)) ./ norm(Z(x2_ind{end}(5:8)));
end

% get indices for states, inputs, and time steps
function [x_ind, u_ind, t_ind, x2_ind, u2_ind, t2_ind] = get_inds(vals)

n = vals.n; m = vals.m; N = vals.N; n_r = vals.n_r; mode = vals.mode;

x_ind = cell(N,1);
u_ind = cell(N-1,1);
t_ind = cell(N-1,1);

x2_ind = cell(N,1);
u2_ind = cell(N-1,1);
t2_ind = cell(N-1,1);

for i = 1:N-1
    x_ind{i} = 1+(i-1)*(n{1}+m+1):n{1}+(i-1)*(n{1}+m+1);
    u_ind{i} = n{1}+1+(i-1)*(n{1}+m+1):n{1}+m+(i-1)*(n{1}+m+1);
    t_ind{i} = n{1}+m+1+(i-1)*(n{1}+m+1);
    
    x2_ind{i} = 1+(i-1)*(n{mode}+m+1):n{mode}+(i-1)*(n{mode}+m+1);
    u2_ind{i} = n{mode}+1+(i-1)*(n{mode}+m+1):n{mode}+m+(i-1)*(n{mode}+m+1);
    t2_ind{i} = n{mode}+m+1+(i-1)*(n{mode}+m+1);
end

x_ind{N} = 1+(N-1)*(n{1}+m+1):n{1}+(N-1)*(n{1}+m+1);
x2_ind{N} = 1+(N-1)*(n{mode}+m+1):n{mode}+(N-1)*(n{mode}+m+1);

end

% discrete dynamics
function xk1 = discrete_dynamics(model, x, u, dt)

tmp = model.discrete_dynamics('x', x, 'u', u, 'dt', dt);
xk1 = full(tmp.xk1);

end

% discrete dynamics Jacobian
function [Ak, Bk, Tk] = discrete_dynamics_jacobian(model, x, u, dt)

tmp = model.discrete_dynamics_jacobian('x', x, 'u', u, 'dt', dt);
Ak = full(tmp.Ak);
Bk = full(tmp.Bk);
Tk = full(tmp.Tk);

end

% operator for quaternion error
function L = lmult(q)

v_hat = [0, -q(4), q(3); q(4), 0, -q(2); -q(3), q(2), 0];
L = [q(1), -q(2:4)'; q(2:4), q(1)*eye(3)+v_hat];

end

% calculate state error with 3-parameter expression for quaternions
function dx = state_error(x,x0,vals)

dx = zeros(vals.n_r,1);
q_diff_b = lmult(x0(1:4))'*x(1:4);
q_diff_t = lmult(x0(5:8))'*x(5:8);

dx(1:3) = 1/q_diff_b(1)*q_diff_b(2:4);
dx(4:6) = 1/q_diff_t(1)*q_diff_t(2:4);

if vals.mode == 1
    dx(7:end) = x(9:14) - x0(9:14);
elseif vals.mode == 2
    dx(7:end) = x(12:17) - x0(12:17);
end

end

% state error Jacobian to reduced quaternion coordinates
function J = state_error_jacobian(x,vals)

n = vals.n; m = vals.m; N = vals.N; n_r = vals.n_r;

J = zeros(n{vals.mode},n_r);

q_b = x(1:4);
v_hat_b = [0, -q_b(4), q_b(3); q_b(4), 0, -q_b(2); -q_b(3), q_b(2), 0];
G_b = [-q_b(2:4)'; q_b(1)*eye(3) + v_hat_b];

q_t = x(5:8);
v_hat_t = [0, -q_t(4), q_t(3); q_t(4), 0, -q_t(2); -q_t(3), q_t(2), 0];
G_t = [-q_t(2:4)'; q_t(1)*eye(3) + v_hat_t];

J(1:4,1:3) = G_b;
J(5:8,4:6) = G_t;

if vals.mode == 1
    J(9:14,7:n_r) = eye(6);
elseif vals.mode == 2
    J(12:17,7:n_r) = eye(6);
end

end