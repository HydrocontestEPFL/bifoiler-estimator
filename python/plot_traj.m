data = csvread('../build/sim_x.csv',1,0);

rx = data(:,7);
ry = data(:,8);
rz = data(:,9);

close all; hold off;
plot3(rx, ry, rz)
hold on;
grid on;
plot3(rx, ry, rz(1)*ones(length(rz),1), ':k')

daspect([1 1 1])
title('Simulation trajectory')
xlabel('r_x [m]')
ylabel('r_y [m]')
zlabel('r_z [m]')
legend('trajectory', 'projection-xy')
printpdf(gcf, 'trajectory_3d_matlab.pdf', 1, 0.7)

