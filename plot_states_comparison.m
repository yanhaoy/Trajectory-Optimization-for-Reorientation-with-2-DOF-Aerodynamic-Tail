function plot_states_comparison(Z_1,Z_2,Z_3,x_ind,u_ind,t_ind,vals)

n = vals.n; m = vals.m; N = vals.N;

t_vec = zeros(N,3);
rb_vec = zeros(N,3);
pb_vec = zeros(N,3);
y_vec = zeros(N,3);
rt_vec = zeros(N,3);
pt_vec = zeros(N,3);
u1_vec = zeros(N-1,3);
u2_vec = zeros(N-1,3);

for i = 1:N
    
    for j = 1:3
        
        if j == 1
            Z = Z_1;
        elseif j == 2
            Z = Z_2;
        else
            Z = Z_3;
        end
        
        x = Z(x_ind{i});

        eul_b = quat2eul(x(1:4)');
        eul_t = quat2eul(x(5:8)');

        rb_vec(i,j) = eul_b(1);
        pb_vec(i,j) = eul_b(2);
        y_vec(i,j) = eul_b(3);
        rt_vec(i,j) = eul_t(1);
        pt_vec(i,j) = eul_t(2);
        
        if i < N
            u = Z(u_ind{i});

            u1_vec(i,j) = u(1);
            u2_vec(i,j) = u(2);
            
            t = Z(t_ind{i});

            t_vec(i+1,j) = t_vec(i,j) + t;
        end
        
    end
    
end

hold on

subplot(2,4,1)
hold on
inert = plot(t_vec(:,1),rb_vec(:,1),'b');
heavy = plot(t_vec(:,2),rb_vec(:,2),'c');
light = plot(t_vec(:,3),rb_vec(:,3),'r');
title('body roll')

subplot(2,4,2)
hold on
plot(t_vec(:,1),pb_vec(:,1),'b');
plot(t_vec(:,2),pb_vec(:,2),'c');
plot(t_vec(:,3),pb_vec(:,3),'r');
title('body pitch')

subplot(2,4,3)
hold on
plot(t_vec(:,1),y_vec(:,1),'b');
plot(t_vec(:,2),y_vec(:,2),'c');
plot(t_vec(:,3),y_vec(:,3),'r');
title('body/tail yaw')

subplot(2,4,4)
hold on
plot(t_vec(:,1),rt_vec(:,1),'b');
plot(t_vec(:,2),rt_vec(:,2),'c');
plot(t_vec(:,3),rt_vec(:,3),'r');
title('tail roll')

subplot(2,4,5)
hold on
plot(t_vec(:,1),pt_vec(:,1),'b');
plot(t_vec(:,2),pt_vec(:,2),'c');
plot(t_vec(:,3),pt_vec(:,3),'r');
title('tail pitch')

subplot(2,4,6)
hold on
plot(t_vec(1:end-1,1),u1_vec(:,1),'b');
plot(t_vec(1:end-1,2),u1_vec(:,2),'c');
plot(t_vec(1:end-1,3),u1_vec(:,3),'r');
title('input 1')

subplot(2,4,7)
hold on
plot(t_vec(1:end-1,1),u2_vec(:,1),'b');
plot(t_vec(1:end-1,2),u2_vec(:,2),'c');
plot(t_vec(1:end-1,3),u2_vec(:,3),'r');
title('input 2')

hL = legend([inert,heavy,light],{'Inertial','Aerodynamic (Heavy)','Aerodynamic (Light)'},'Location','best');