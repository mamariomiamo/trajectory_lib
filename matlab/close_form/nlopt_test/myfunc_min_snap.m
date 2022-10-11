function [val] = myfunc_min_snap(x, waypts, n_order, r_order, k_segment, coefficient_number, dim, time_penalty)
% free decision variable
free_number = r_order*(k_segment-1); % intermedia v,a
d_p = zeros(free_number, dim); % vector to store free decision variables

opt_var_num = free_number*dim + k_segment;

time_vector = x((opt_var_num-k_segment+1):end);
time_vector = [0,time_vector];
derivative_vector = x(1:opt_var_num-k_segment);
d_p = reshape(derivative_vector, [free_number,dim]);

mapping_A=getMappingA(n_order, k_segment, time_vector, coefficient_number);
fixed_number = 2*(r_order + 1) + k_segment-1; % initial and final p,v,a + intermediate p
d_f = zeros(fixed_number, dim); % vector to store fixed decision variables


C = getSelectionC(fixed_number, free_number, k_segment, r_order);
C_t = C';
Q_hessian = calQ(time_vector);
M = mapping_A;
d_f(1,:) =  waypts(:,1);
for i=4:4+k_segment-1
    d_f(i,:) = waypts(:,i-2);
end
R = C * inv(M)' * Q_hessian * inv(M) * C_t;
d = [d_f;d_p];
traj_cost = d(:,1)' * R * d(:,1) + d(:,2)' * R * d(:,2) + d(:,3)' * R * d(:,3);
time_cost = time_vector(end) * time_penalty;
val = traj_cost + time_cost;