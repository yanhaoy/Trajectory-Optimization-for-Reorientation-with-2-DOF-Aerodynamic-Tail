clear; clc; close all

% reorienting body with inertial tail
%addpath('utils')
%addpath('utils\casadi')
%import casadi.*

% state and input sizes
vals.n = 11;
vals.m = 5;

% set time span and knot points
vals.t0 = 0;
vals.tf = 10;
vals.N = 100;
vals.dt = (vals.tf-vals.t0) / vals.N;
vals.Nh = vals.n*vals.N + vals.m*(vals.N-1);
[x_ind, u_ind] = get_inds(vals);

% initial and final conditions
vals.x0 = zeros(vals.Nh,1); 
R0 = pi*rand(3,1) - pi/2;
for i = 1:vals.N
    vals.x0(1+(i-1)*(vals.n+vals.m):3+(i-1)*(vals.n+vals.m)) = R0;
end
vals.xf = zeros(vals.Nh,1);

% using Euler angles
addpath('euler_angle_2dof')

model = init_dynamics();

% weight matrices
P = eye(vals.n);
R = 2*eye(vals.m);
H = diag([repmat([diag(P)',diag(R)'],1,vals.N-1),diag(P)']);

% cost function
cost_fun = @(x)x'*H*x;

% solve constrained optimization
options = optimoptions('fmincon','SpecifyConstraintGradient',true,'Display','iter','MaxFunctionEvaluations',1e4);
Z = fmincon(cost_fun,vals.x0,[],[],[],[],[],[],@(x)nonlcon(x,model,vals),options);

figure();
viz_tail(Z,x_ind,u_ind,vals);
figure();
plot_states(Z,x_ind,u_ind,vals)

% constaints
function [cineq,ceq,Gineq,Geq]=nonlcon(Z,model,vals)

n = vals.n; m = vals.m; N = vals.N;
[x_ind, u_ind] = get_inds(vals);

ceq = zeros(n*(N+1)+3*(N-1),1);
Geq = zeros(vals.Nh,n*(N+1));
cineq = [];
Gineq = [];

for i = 1:N-1
    x_dyn = discrete_dynamics(model,...
        Z(x_ind{i}),Z(u_ind{i}),vals.dt);
    ceq(n*(i-1)+1:n*(i-1)+n) = x_dyn - Z(x_ind{i+1});
    
    [Ak, Bk] = discrete_dynamics_jacobian(model, Z(x_ind{i}), Z(u_ind{i}), vals.dt);
    
    Geq(x_ind{i},n*(i-1)+1:n*(i-1)+n) = Ak';
    Geq(u_ind{i},n*(i-1)+1:n*(i-1)+n) = Bk';
    Geq(x_ind{i+1},n*(i-1)+1:n*(i-1)+n) = -eye(n);
    
    ceq(n*(N+1)+3*(i-1)+1:n*(N+1)+3*(i-1)+3) = Z(u_ind{i}(1:3));
    Geq(u_ind{i}(1:3),n*(N+1)+3*(i-1)+1:n*(N+1)+3*(i-1)+3) = eye(3);
end

ceq(n*(N-1)+1:n*(N-1)+n) = Z(x_ind{1}) - vals.x0(1:n);
Geq(x_ind{1},n*(N-1)+1:n*(N-1)+n) = eye(n);

ceq(n*N+1:n*N+n) = Z(x_ind{end}) - vals.xf(1:n);
Geq(x_ind{end},n*N+1:n*N+n) = eye(n);
end

function [x_ind, u_ind] = get_inds(vals)

n = vals.n; m = vals.m; N = vals.N;

x_ind = cell(N,1);
u_ind = cell(N,1);

for i = 1:N-1
    x_ind{i} = 1+(i-1)*(n+m):n+(i-1)*(n+m);
    u_ind{i} = n+1+(i-1)*(n+m):n+m+(i-1)*(n+m);
end
x_ind{N} = 1+(N-1)*(n+m):n+(N-1)*(n+m);
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