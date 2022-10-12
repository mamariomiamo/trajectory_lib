function [derivative, C, d_f, fixed_number, free_number] = minSnapCloseFormOpt(waypts, time_vector, dim, r_order, n_order, k_segment, coefficient_number)
% % A is the mapping matrix between p (coefficients) and d (derivatives)
mapping_A=getMapping(n_order, k_segment, time_vector, coefficient_number);

% fixed decision variable
fixed_number = 2*(r_order + 1) + k_segment-1; % initial and final p,v,a + intermediate p
d_f = zeros(fixed_number, dim); % vector to store fixed decision variables

% free decision variable
free_number = r_order*(k_segment-1); % intermedia v,a
C = getSelectionC(fixed_number, free_number, k_segment, r_order);
C_t = C';

% time_vector = allocate_time(waypts, T);
% time_vector = [0 2 4 10];
Q_hessian = getHessian(time_vector);
M = mapping_A;
R = C * inv(M)' * Q_hessian * inv(M) * C_t;
% R = C * inv(M)' * Q_hessian) \M * C_t;

[R_ff, R_fp, R_pf, R_pp]=getSubR(R, fixed_number, free_number);

% d_primal_x = getDerivativePrimal(d_f(:,1), R_pp, R_fp, k_segment, waypts(1,:));
% d_primal_y = getDerivativePrimal(d_f(:,2), R_pp, R_fp, k_segment, waypts(2,:));
% d_primal_z = getDerivativePrimal(d_f(:,3), R_pp, R_fp, k_segment, waypts(3,:));

d_primal_initial_guess = getDerivativePrimal_3axis(d_f, R_pp, R_fp, k_segment, waypts);

derivative = d_primal_initial_guess(fixed_number+1:end,:);
end