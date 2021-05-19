clear; clc; close all

addpath('utils')
addpath('utils\casadi')
import casadi.*

% reference state and input sizes
vals.n{1} = 14; % full state space
vals.n{2} = 20; % full state space with body COM
vals.m = 4; % input space with 2 slack variables

vals.mode = 1;

% set time span and knot points
vals.t0 = 0;
vals.tf = 5;
vals.N = 100;
vals.dt = (vals.tf-vals.t0) / vals.N; % nominal time step

% number of primals
vals.Nh{1} = vals.n{1}*vals.N + vals.m*(vals.N-1) + vals.N-1;
vals.Nh{2} = vals.n{2}*vals.N + vals.m*(vals.N-1) + vals.N-1;

% get indices from packed states
% first set mainly for plotting; second set used in optimization
[x_ind, u_ind, t_ind, ~, ~, ~] = get_inds(vals);

% reference trajectory
vals.xref{1} = zeros(vals.Nh{1},1); 

vals.R0 = [pi/6 pi/6 pi/6]; %pi/2*rand(1,3) - pi/4; % random initial orientation
vals.Rf = [0,0,0]; % final orientation
for i = 1:vals.N
    vals.xref{1}(x_ind{i}(1:4)) = eul2quat(vals.Rf,'XYZ')';
    vals.xref{1}(x_ind{i}(5:8)) = eul2quat(vals.Rf,'XYZ')';
    
    if i < vals.N
        vals.xref{1}(t_ind{i}) = vals.dt;
        vals.xref{1}(u_ind{i}(3:4)) = 1;
    end
end

%% using inertial tail with quaternions
addpath('quaternion')
model_inert = init_dynamics();
vals.mode = 1;

Z = SX.sym('Z', size(vals.xref{vals.mode}));
f = cost_fun(Z, vals);
[g, lbg, ubg] = nonlcon(Z,model_inert,vals);
[lbx, ubx] = get_bounds(vals);

% Create an NLP solver
prob = struct('f', f, 'x', Z, 'g', g);
opts = struct;
opts.ipopt.max_cpu_time = 1e4;
opts.ipopt.max_iter = 1e6;
opts.ipopt.hessian_approximation = 'limited-memory';
solver = nlpsol('solver', 'ipopt', prob, opts);

% Solve the NLP
% load Z_inert.mat
Z_inert = vals.xref{vals.mode};
sol = solver('x0', Z_inert, 'lbx', lbx, 'ubx', ubx,...
            'lbg', lbg, 'ubg', ubg);
Z_inert = full(sol.x);

save Z_inert.mat Z_inert

%% using heavy aerodynamic tail
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
        vals.xref{2}(u2_ind{i}) = vals.xref{1}(u_ind{i});
    end
end

%%
Z = casadi.SX.sym('Z', size(vals.xref{vals.mode}));
f = cost_fun(Z, vals);
[g, lbg, ubg] = nonlcon(Z,model_heavy,vals);
[lbx, ubx] = get_bounds(vals);

% Create an NLP solver
prob = struct('f', f, 'x', Z, 'g', g);
opts = struct;
% opts.ipopt.max_cpu_time = 1800;
% opts.ipopt.max_iter = 1e6;
opts.ipopt.hessian_approximation = 'limited-memory';
solver = casadi.nlpsol('solver', 'ipopt', prob, opts);

% Solve the NLP
% load Z_heavy.mat
Z_heavy = vals.xref{vals.mode};
sol = solver('x0', Z_heavy, 'lbx', lbx, 'ubx', ubx,...
            'lbg', lbg, 'ubg', ubg);
Z_heavy = full(sol.x);

save Z_heavy.mat Z_heavy

%% using light aerodynamic tail
addpath('quaternion_aerodynamics_free_falling')
model_light = init_dynamics_light();
vals.mode = 2;

Z = casadi.SX.sym('Z', size(vals.xref{vals.mode}));
f = cost_fun(Z, vals);
[g, lbg, ubg] = nonlcon(Z,model_light,vals);
[lbx, ubx] = get_bounds(vals);

% Create an NLP solver
prob = struct('f', f, 'x', Z, 'g', g);
opts = struct;
% opts.ipopt.max_cpu_time = 1800;
% opts.ipopt.max_iter = 1e6;
opts.ipopt.hessian_approximation = 'limited-memory';
solver = casadi.nlpsol('solver', 'ipopt', prob, opts);

% Solve the NLP
% load Z_light.mat
Z_light = vals.xref{vals.mode};
sol = solver('x0', Z_light, 'lbx', lbx, 'ubx', ubx,...
            'lbg', lbg, 'ubg', ubg);
Z_light = full(sol.x);

save Z_light.mat Z_light

%%
%[cineq,ceq,Gineq,Geq]=nonlcon(Z,model,vals);
%viz_tail(Z,x_ind,u_ind,vals);

%%
figure('Name','Intertial Tail States');
plot_states_comparison(Z_inert,Z_heavy,Z_light,x_ind,u_ind,t_ind,x2_ind,u2_ind,t2_ind,vals);

% cost function
function [f] = cost_fun(Z, vals)

N = vals.N;

[~,~,~, ~, ~, t2_ind] = get_inds(vals);

f = 0;

for i = 1:N-1
   f = f + Z(t2_ind{i}); 
end

end

% constraints
function [g, lbg, ubg] = nonlcon(Z,model,vals)

n = vals.n{vals.mode}; N = vals.N;
[~,~,~, x2_ind, u2_ind, t2_ind] = get_inds(vals);

ceq = [];

for i = 1:N-1
    x = Z(x2_ind{i}); 
    u = Z(u2_ind{i}(1:2));
    t = Z(t2_ind{i});
    s = Z(u2_ind{i}(3:4));
    
    x_dyn = discrete_dynamics(model,x,u,t);
    x_norm = [x_dyn(1:4)*s(1); x_dyn(5:8)*s(2); x_dyn(9:end)];
    
    % dynamics contraints
    ceq = [ceq; state_error(x_norm,Z(x2_ind{i+1}),vals)];
    
    % unit quaternion constraints
    ceq = [ceq; x(1:4)'*x(1:4) - 1; x(5:8)'*x(5:8) - 1];
    
end

% final unit quaternion constraint
x = Z(x2_ind{end}); 
ceq = [ceq; x(1:4)'*x(1:4) - 1; x(5:8)'*x(5:8) - 1];

% initial state constraint
x_int = [eul2quat(vals.R0,'XYZ')'; eul2quat(vals.Rf,'XYZ')'; zeros(n-8, 1)];
% ceq = [ceq; Z(x2_ind{1}) - x_int];
ceq = [ceq; state_error(Z(x2_ind{1}),x_int,vals)];
% ceq = [ceq; 1 - abs(Z(x2_ind{1}(1:4))'*x_int(1:4)); 1 - abs(Z(x2_ind{1}(5:8))'*x_int(5:8))];

% final state constraint
x_term = [eul2quat(vals.Rf,'XYZ')'; eul2quat(vals.Rf,'XYZ')'; zeros(n-8, 1)];
% ceq = [ceq; Z(x2_ind{end}) - x_term];
g_term = state_error(Z(x2_ind{end}),x_term,vals);
if vals.mode == 2 
    ceq = [ceq; g_term([1:6, 10:15])];
else
    ceq = [ceq; g_term];
end
% ceq = [ceq; 1 - abs(Z(x2_ind{end}(1:4))'*x_term(1:4)); 1 - abs(Z(x2_ind{end}(5:8))'*x_term(5:8))];

g = ceq;
lbg = zeros(size(ceq));
ubg = zeros(size(ceq));

end

function [lbx, ubx] = get_bounds(vals)

N = vals.N;
[~,~,~, x2_ind, u2_ind, t2_ind] = get_inds(vals);

lbx = -inf*ones(vals.Nh{vals.mode}, 1);
ubx = inf*ones(vals.Nh{vals.mode}, 1);

for i = 1:N-1
    lbx(x2_ind{i}(1:8)) = -1;
    ubx(x2_ind{i}(1:8)) = 1;
    
    lbx(u2_ind{i}(1:2)) = -1;
    ubx(u2_ind{i}(1:2)) = 1;
    
    lbx(u2_ind{i}(3:4)) = 0.9;
    ubx(u2_ind{i}(3:4)) = 1.1;
    
    lbx(t2_ind{i}) = 0.001;
    ubx(t2_ind{i}) = 0.05;
end

lbx(x2_ind{end}(1:8)) = -1;
ubx(x2_ind{end}(1:8)) = 1;

% initial state constraint
% if vals.mode == 1 
%     x_int = [eul2quat(vals.R0,'XYZ')'; eul2quat(vals.Rf,'XYZ')'; zeros(6, 1)];
% else
%     x_int = [eul2quat(vals.R0,'XYZ')'; eul2quat(vals.Rf,'XYZ')'; zeros(12, 1)];
% end
% lbx(x2_ind{1}(9:end)) = x_int(9:end);
% ubx(x2_ind{1}(9:end)) = x_int(9:end);

% final state constraint
% if vals.mode == 1 
%     x_term = [eul2quat(vals.Rf,'XYZ')'; eul2quat(vals.Rf,'XYZ')'; zeros(6, 1)];
% else
%     x_term = [eul2quat(vals.Rf,'XYZ')'; eul2quat(vals.Rf,'XYZ')'; zeros(12, 1)];
% end
% lbx(x2_ind{end}(9:end)) = x_term(9:end);
% ubx(x2_ind{end}(9:end)) = x_term(9:end);

end

% get indices for states, inputs, and time steps
function [x_ind, u_ind, t_ind, x2_ind, u2_ind, t2_ind] = get_inds(vals)

n = vals.n; m = vals.m; N = vals.N; mode = vals.mode;

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
xk1 = tmp.xk1;

end

% operator for quaternion error
function L = lmult(q)

v_hat = [0, -q(4), q(3); q(4), 0, -q(2); -q(3), q(2), 0];
L = [q(1), -q(2:4)'; q(2:4), q(1)*eye(3)+v_hat];

end

% calculate state error with 3-parameter expression for quaternions
function dx = state_error(x,x0,vals)

q_diff_b = lmult(x0(1:4))'*x(1:4);
q_diff_t = lmult(x0(5:8))'*x(5:8);

dx = [1/q_diff_b(1)*q_diff_b(2:4); 1/q_diff_t(1)*q_diff_t(2:4); x(9:end) - x0(9:end)];

end
