clear, clc, close all;
r_order = 2; %p,v,a
% r_order = 3; %p,v,a,j

n_order = 5;
coefficient_number = n_order+1;

% A_sub = zeros(dim,dim);

% waypts = [0,0;
%     1,2;
%     2,-1;
%     4,8;
%     5,2]';

waypts = [ -2 4 1;
    -1 0.5 2;
    0 0 3;
    1 1 1]';

k_segment = size(waypts,2) - 1; % 4
% A = zeros(dim*k_segment, dim*k_segment);

T = 1.78179;
% T = 10;
time_vector = allocate_time(waypts,T);
% time_vector = [0 2 4 10];
% A is the mapping matrix between p (coefficients) and d (derivatives)
mapping_A=getMappingA(n_order, k_segment, time_vector, coefficient_number);

% fixed decision variable
fixed_number = 2*(r_order + 1) + k_segment-1; % initial and final p,v,a + intermediate p
d_f = zeros(fixed_number, 1); % vector to store fixed decision variables

% free decision variable
free_number = r_order*(k_segment-1); % intermedia v,a
d_p = zeros(free_number, 1); % vector to store free decision variables

% C_t * [d_f; d_p] = d
C = getSelectionC(fixed_number, free_number, k_segment, r_order);
C_t = C';

% time_vector = allocate_time(waypts, T);
% time_vector = [0 2 4 10];
Q_hessian = calQ(time_vector);
M = mapping_A;
R = C * inv(M)' * Q_hessian * inv(M) * C_t;
% R = C * inv(M)' * Q_hessian) \M * C_t;

[R_ff, R_fp, R_pf, R_pp]=getSubR(R, fixed_number, free_number);

d_primal_x = getDerivativePrimal(d_f, R_pp, R_fp, k_segment, waypts(1,:));
d_primal_y = getDerivativePrimal(d_f, R_pp, R_fp, k_segment, waypts(2,:));
d_primal_z = getDerivativePrimal(d_f, R_pp, R_fp, k_segment, waypts(3,:));

% map d_primal to full set of derivatives using the selection matrix
%
p_x = inv(M)  * C_t * d_primal_x;
p_y = inv(M)  * C_t * d_primal_y;
p_z = inv(M)  * C_t * d_primal_z;

p_x = reshape(p_x,[n_order+1,k_segment]);
p_y = reshape(p_y,[n_order+1,k_segment]);
p_z = reshape(p_z,[n_order+1,k_segment]);

[f1, x, vx, ax] = plotTrajectory1D(n_order, k_segment, time_vector, p_x, waypts(1,:), 1, 'x axis');
movegui(f1, 'northwest');
[f2, y, vy, ay] = plotTrajectory1D(n_order, k_segment, time_vector, p_y, waypts(2,:), 2, 'y axis');
movegui(f2, 'north');
[f3, z, vz, az] = plotTrajectory1D(n_order, k_segment, time_vector, p_z, waypts(3,:), 3, 'z axis');
movegui(f3, 'northeast');

f4 = plotTrajectory2D(x,y,waypts,4);
movegui(f4, 'southeast');

f5 = plotTrajectory3D(x,y,z,waypts,5);
movegui(f5, 'southwest');