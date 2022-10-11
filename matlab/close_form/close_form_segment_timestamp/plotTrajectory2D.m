function f = plotTrajectory2D(x,y,waypts,fig_index)
f = figure(fig_index);
hold on;
plot(waypts(1,:), waypts(2,:),'b--');
plot(waypts(1,:), waypts(2,:),'*r');
plot(x,y,'LineWidth',3);
end