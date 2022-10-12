function [f, x, vx, ax, t] = plotTrajectory1axis(n_order, k_segment, time_vector, p_x, waypts, fig_index, fig_name)
f = figure(fig_index);
sgtitle(fig_name);
subplot(3,1,3);
subplot(3,1,1);
time = plot_time(time_vector);
plot(time,waypts,'*r');
hold on;
% subplot(3,1,1);
% plot(time_vector,waypts(1,:),'b--');
% each column of p is the coefficient vector for one segment
x = [];
vx = [];
ax = [];
y = [];
t = [];
for i=1:k_segment
    for t_tick=0:0.01:time_vector(i)
        x = [x poly_evaluate(0, t_tick, n_order) * p_x(:,i)];
        vx = [vx poly_evaluate(1, t_tick, n_order) * p_x(:,i)];
        ax = [ax poly_evaluate(2, t_tick, n_order) * p_x(:,i)];
        t = [t t_tick+time(i)];
    end
end

% t_segment=time_vector(1):0.01:time_vector(k_segment+1);
subplot(3,1,1);
plot(t,x,'r')
title("Position")
% figure(2)
% hold on;
for i=1:k_segment
    xline(time(i),'k--');
end
subplot(3,1,2);
plot(t,vx,'r')
title("Velocity")
% figure(3)
% hold on;
for i=1:k_segment
    xline(time(i),'k--');
end
subplot(3,1,3);
plot(t,ax,'r')
title("Acceleration")
for i=1:k_segment
    xline(time(i),'k--');
end
end

