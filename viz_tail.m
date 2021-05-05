function viz_tail(Z,x_ind,u_ind,vals)

n = vals.n; m = vals.m; N = vals.N;

bRecord = 1;  % Uncomment this to save a video
if bRecord
    % Define video recording parameters
    v = VideoWriter('tail_reorientation', 'MPEG-4');
    v.Quality = 100;
    v.FrameRate = 30;
    open(v);
end

for i = 1:N
    clf
    hold on
    axis([-1.5,1.5,-1.5,1.5,-1.5,1.5])
    view(-45,-60)
    
    x = Z(x_ind{i});
    R_b = quat2rotm(x(1:4)');
    R_t = quat2rotm(x(5:8)')*quat2rotm(x(1:4)');
    
    scatter3(0,0,0,1500,'k','filled');
    plot3([0,R_b(1,1)],[0,R_b(2,1)],[0,R_b(3,1)],'k');
    plot3([0,R_b(1,2)],[0,R_b(2,2)],[0,R_b(3,2)],'r');
    plot3([0,R_b(1,3)],[0,R_b(2,3)],[0,R_b(3,3)],'g');
    
    scatter3(R_b(1,1),R_b(2,1),R_b(3,1),150,[0.3,0.3,0.3],'filled');
    plot3([R_b(1,1),0.2*R_t(1,1)+R_b(1,1)],[R_b(2,1),0.2*R_t(2,1)+R_b(2,1)],[R_b(3,1),0.2*R_t(3,1)+R_b(3,1)],'b');
    plot3([R_b(1,1),0.2*R_t(1,2)+R_b(1,1)],[R_b(2,1),0.2*R_t(2,2)+R_b(2,1)],[R_b(3,1),0.2*R_t(3,2)+R_b(3,1)],'r');
    plot3([R_b(1,1),0.2*R_t(1,3)+R_b(1,1)],[R_b(2,1),0.2*R_t(2,3)+R_b(2,1)],[R_b(3,1),0.2*R_t(3,3)+R_b(3,1)],'g');
    
     if bRecord
        frame = getframe(gcf);
        writeVideo(v,frame);
     end
end

if bRecord
    close(v);
end