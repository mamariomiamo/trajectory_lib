function f = plotTrajectory3D(x,y,z,waypts,fig_index)
f = figure(fig_index);
plot3(waypts(1,:), waypts(2,:), waypts(3,:),'b--');
hold on
% grid on
plot3(waypts(1,:), waypts(2,:), waypts(3,:),'*r');
plot3(x,y,z,'LineWidth',3);
xlabel('x position') 
ylabel('y position')
zlabel('z position') 
end