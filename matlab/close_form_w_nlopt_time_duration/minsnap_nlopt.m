% min snap we dun have nonlinear constraints
% only needs lb and ub, f, f_data, xinit
% call [xopt, fmin, retcode] = nlopt_minimize(algorithm, f, f_data, lb, ub,
%                                               xinit, stop)
clc, clear, close all;
r_order = 2; %p,v,a
n_order = 5;
coefficient_number = n_order+1;

dim = 3; % 3 as in x,y,z
v_max_abs = 5.0;
a_max_abs = 5.0;
v_nominal = v_max_abs;
kinematic_constraint_low = [-5 -5];
kinematic_constraint_up = [5, 5];
% waypts = [0,0, 0;
%     1,2, 1;
%     2,-1, 2;
%     4,8, 3;
%     5,2, 0]';

waypts = [ -20 4 1;
           -10 0.5 2;
            50 30 3;
            100 10 1]';

k_segment = size(waypts,2) - 1;

total_dist = calDistance(waypts, k_segment);

T = total_dist/1;

% time_initial_guess = allocate_time(waypts,T);
time_initial_guess = [0.1 0.1 0.1];
% C for selection matrix
[derivative_initial_guess, C, d_f, fixed_number, free_number] = minSnapCloseFormOpt(waypts, time_initial_guess, dim, r_order, n_order, k_segment, coefficient_number);
derivative_initial_guess = reshape(derivative_initial_guess, [1, free_number*dim]);
derivative_initial_guess = checkFeasibility(derivative_initial_guess,v_max_abs, a_max_abs,k_segment);
derivative_initial_guess = ones(1,length(derivative_initial_guess))*1;
initial_guess = [derivative_initial_guess, time_initial_guess];
C_t = C';

%% optimization lower and upper bounds
% optimizing_variable 
% number = (k_segment-1)*r_order*dim + k_segment
% x = [dvx1,dax1,... dvx(M-1),dax(M-1),dvy1,dyx1,... dvy(M-1),day(M-1), dvz1,dyz1,... dvz(M-1),daz(M-1), t1, ..., tM]
opt_var_num = (k_segment-1)*r_order*dim + k_segment;
derivative_lower_bounds = ones((k_segment-1)*dim,r_order) .* kinematic_constraint_low;
derivative_lower_bounds = reshape(derivative_lower_bounds,[1,(k_segment-1)*r_order*dim])

derivative_upper_bounds = ones((k_segment-1)*dim,r_order) .* kinematic_constraint_up;
derivative_upper_bounds = reshape(derivative_upper_bounds,[1,(k_segment-1)*r_order*dim])

time_lower_bounds = zeros(1,k_segment);
time_upper_bounds = Inf(1,k_segment);

lower_bounds = [derivative_lower_bounds, time_lower_bounds];
upper_bounds = [derivative_upper_bounds, time_upper_bounds];

time_penalty = 10;

%% Setup optimizer
opt.algorithm = NLOPT_LN_COBYLA;
opt.verbose = 1;
opt.lower_bounds = lower_bounds;
opt.upper_bounds = upper_bounds;
opt.min_objective = @(x) myfunc_min_snap(x, waypts, n_order, r_order, k_segment, coefficient_number, dim, time_penalty);
opt.fc = {};

opt.xtol_rel = 0.001;

[xopt, fmin, retcode] = nlopt_optimize(opt, initial_guess);

d_primal = xopt(1:length(xopt)-k_segment);
time_primal = xopt(length(xopt)-(k_segment-1):end);

d_primal = reshape(d_primal,[r_order*(k_segment-1),dim]);
d = [d_f;d_primal];

mapping_matrix_primal=getMapping(n_order, k_segment, time_primal, coefficient_number);
p_primal = inv(mapping_matrix_primal)  * C_t * d;

p_primal = reshape(p_primal, [n_order+1, k_segment,dim]);

total_duration = time_primal * ones(k_segment,1)

p_x = p_primal(:,:,1);
p_y = p_primal(:,:,2);
p_z = p_primal(:,:,3);

[f1, x, vx, ax, tx] = plotTrajectory1axis(n_order, k_segment, time_primal, p_x, waypts(1,:), 1, 'x axis');
movegui(f1, 'northwest');
[f2, y, vy, ay, ty] = plotTrajectory1axis(n_order, k_segment, time_primal, p_y, waypts(2,:), 2, 'y axis');
movegui(f2, 'north');
[f3, z, vz, az, tz] = plotTrajectory1axis(n_order, k_segment, time_primal, p_z, waypts(3,:), 3, 'z axis');
movegui(f3, 'northeast');

f4 = plotTrajectory2axis(x,y,waypts,4);
movegui(f4, 'southeast');

f5 = plotTrajectory3axis(x,y,z,waypts,5);
movegui(f5, 'southwest');