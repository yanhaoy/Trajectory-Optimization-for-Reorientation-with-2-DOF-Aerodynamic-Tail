function plot_states(Z,x_ind,u_ind,t_ind,vals)

n = vals.n; m = vals.m; N = vals.N;

t_vec = zeros(N,1);
rb_vec = zeros(N,1);
pb_vec = zeros(N,1);
y_vec = zeros(N,1);
rt_vec = zeros(N,1);
pt_vec = zeros(N,1);
u1_vec = zeros(N-1,1);
u2_vec = zeros(N-1,1);

for i = 1:N
    x = Z(x_ind{i});
    if n == 11
        rb_vec(i) = x(1);
        pb_vec(i) = x(2);
        y_vec(i) = x(3);
        rt_vec(i) = x(4);
        pt_vec(i) = x(5);
        if i < N
            u = Z(u_ind{i});
            u1_vec(i) = u(4);
            u2_vec(i) = u(5);
            
            t = Z(t_ind{i});
            t_vec(i+1) = t_vec(i) + t;
        end
    else
        eul_b = quat2eul(x(1:4)');
        eul_t = quat2eul(x(5:8)');
        
        rb_vec(i) = eul_b(1);
        pb_vec(i) = eul_b(2);
        y_vec(i) = eul_b(3);
        rt_vec(i) = eul_t(1);
        pt_vec(i) = eul_t(2);
        if i < N
            u = Z(u_ind{i});
            u1_vec(i) = u(1);
            u2_vec(i) = u(2);
            
            t = Z(t_ind{i});
            t_vec(i+1) = t_vec(i) + t;
        end
    end
end

subplot(2,4,1)
plot(t_vec,rb_vec);
title('body roll')

subplot(2,4,2)
plot(t_vec,pb_vec);
title('body pitch')

subplot(2,4,3)
plot(t_vec,y_vec);
title('yaw')

subplot(2,4,4)
plot(t_vec,rt_vec);
title('tail roll')

subplot(2,4,5)
plot(t_vec,pt_vec);
title('tail pitch')

subplot(2,4,6)
plot(t_vec(1:end-1),u1_vec);
title('input 1')

subplot(2,4,7)
plot(t_vec(1:end-1),u2_vec);
title('input 2')