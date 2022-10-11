% min sqrt(x2)
% s.t. x2≥0 , x2≥(a1 x1+b1)^3, and x2≥(a2 x1+b2)^3
% a1=2, b1=0, a2=-1, b2=1.
% min snap we dun have nonlinear constraints
% only needs lb and ub, f, f_data, xinit
% call [xopt, fmin, retcode] = nlopt_minimize(algorithm, f, f_data, lb, ub,
%                                               xinit, stop)
clc, clear, close all;

% opt.algorithm = NLOPT_LD_MMA;
r_order = 2; %p,v,a
% r_order = 3; %p,v,a,j

n_order = 5;
coefficient_number = n_order+1;
% A_sub = zeros(dim,dim);
dim = 3; % 3 as in x,y,z
v_max_abs = 5.0;
a_max_abs = 5.0;
v_nominal = 1.5;
kinematic_constraint_low = [-5 -5];
kinematic_constraint_up = [5, 5];
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

total_dist = calDistance(waypts, k_segment);

T = total_dist/v_nominal;

time_initial_guess = allocate_time(waypts,T);
time_vector = time_initial_guess;
time_initial_guess = time_initial_guess(2:end);

% % A is the mapping matrix between p (coefficients) and d (derivatives)
mapping_A=getMappingA(n_order, k_segment, time_vector, coefficient_number);

% fixed decision variable
fixed_number = 2*(r_order + 1) + k_segment-1; % initial and final p,v,a + intermediate p
d_f = zeros(fixed_number, dim); % vector to store fixed decision variables

% free decision variable
free_number = r_order*(k_segment-1); % intermedia v,a
d_p = zeros(free_number, dim); % vector to store free decision variables

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

d_primal_x = getDerivativePrimal(d_f(:,1), R_pp, R_fp, k_segment, waypts(1,:));
d_primal_y = getDerivativePrimal(d_f(:,2), R_pp, R_fp, k_segment, waypts(2,:));
d_primal_z = getDerivativePrimal(d_f(:,3), R_pp, R_fp, k_segment, waypts(3,:));

d_primal = getDerivativePrimal_3axis(d_f, R_pp, R_fp, k_segment, waypts);

derivative_initial_guess = d_primal(fixed_number+1:end,:);
derivative_initial_guess = reshape(derivative_initial_guess, [1, free_number*dim]);

initial_guess = [derivative_initial_guess, time_initial_guess];
% p = inv(M)  * C_t * d_primal;
% 
% % map d_primal to full set of derivatives using the selection matrix
% %
% p_x = inv(M)  * C_t * d_primal_x;
% p_y = inv(M)  * C_t * d_primal_y;
% p_z = inv(M)  * C_t * d_primal_z;
% 
% p_x = reshape(p_x,[n_order+1,k_segment]);
% p_y = reshape(p_y,[n_order+1,k_segment]);
% p_z = reshape(p_z,[n_order+1,k_segment]);
% 
% p = reshape(p, [n_order+1, k_segment,dim]);

%% optimization

% opt.algorithm = NLOPT_LN_BOBYQA;
opt.algorithm = NLOPT_LN_COBYLA;

opt.verbose = 1;
% optimizing_variable 
% number = k_segment*r_order*dim + k_segment
% x = [dvx1,dax1,dvy1,day1,dvz1,daz1, ..., dvxM,daxM,dvyM,dayM,dvzM,dazM, t1, ..., tM]
% x = [dvx1,dax1,... dvxM,daxM,dvy1,dyx1,... dvyM,dayM, dvz1,dyz1,... dvzM,dazM, t1, ..., tM]
opt_var_num = (k_segment-1)*r_order*dim + k_segment;
derivative_lower_bounds = ones((k_segment-1)*dim,r_order) .* kinematic_constraint_low;
derivative_lower_bounds = reshape(derivative_lower_bounds,[1,(k_segment-1)*r_order*dim])

derivative_upper_bounds = ones((k_segment-1)*dim,r_order) .* kinematic_constraint_up;
derivative_upper_bounds = reshape(derivative_upper_bounds,[1,(k_segment-1)*r_order*dim])

time_lower_bounds = zeros(1,k_segment);
time_upper_bounds = Inf(1,k_segment);

lower_bounds = [derivative_lower_bounds, time_lower_bounds];
upper_bounds = [derivative_upper_bounds, time_upper_bounds];

% opt.lower_bounds = [-inf, 0];
% opt.upper_bounds = [inf, inf];
time_penalty = 10;

opt.lower_bounds = lower_bounds;
opt.upper_bounds = upper_bounds;

opt.min_objective = @(x) myfunc_min_snap(x, waypts, n_order, r_order, k_segment, coefficient_number, dim, time_penalty);

% opt.fc = { (@(x) myconstraint_min_snap(x,2,0)), (@(x) myconstraint_min_snap(x,-1,1)) }; % evaluate constraint
opt.fc = {};
% opt.fc_tol = [1e-8, 1e-8];

% opt.xtol_rel = 1e-4;

[xopt, fmin, retcode] = nlopt_optimize(opt, initial_guess);

% [xopt, fmin, retcode] = nlopt_optimize(opt, [0.3 0.3])