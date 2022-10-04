% function A_sub=getAsub(n_order, k, t_alloc, Q_total, waypts)
% A_sub = zeros(n_order,n_order);
% 
% 
% 
% end
r_order = 2; %p,v,a
% r_order = 3; %p,v,a,j

n_order = 5;
dim = n_order+1;
% A_sub = zeros(dim,dim);

% waypts = [0,0;
%     1,2;
%     2,-1;
%     4,8;
%     5,2]';

waypts = [ -2 4;
    -1 0.5;
    0 0;
    1 1;]';

% waypts = [0,0;
%     1,2;
%     2,-1]';

k_segment = size(waypts,2) - 1; % 4
% A = zeros(dim*k_segment, dim*k_segment);

T = 1.78179;
t_alloc = allocate_time(waypts,T);
% 
% A_sub(1,:) = poly_evaluate(0,t_alloc(1),n_order);
% A_sub(2,:) = poly_evaluate(1,t_alloc(1),n_order);
% A_sub(3,:) = poly_evaluate(2,t_alloc(1),n_order);
% A_sub(4,:) = poly_evaluate(0,t_alloc(2),n_order);
% A_sub(5,:) = poly_evaluate(1,t_alloc(2),n_order);
% A_sub(6,:) = poly_evaluate(2,t_alloc(2),n_order);
A=[];

% every segment will have two sets of derivative constriants
% p,v,a for start of the segment
% p,v,a for end of the segment
for i=1:k_segment
    A_sub = zeros(dim);
    for j=1:3
        A_sub(j,:) = poly_evaluate(j-1,t_alloc(i),n_order);
    end
    for j=4:6
        A_sub(j,:) = poly_evaluate(j-4,t_alloc(i+1),n_order);
    end
    A = blkdiag(A, A_sub);
end

% fixed decision variable
fixed_number = 2*(r_order + 1) + k_segment-1; % initial and final p,v,a + intermediate p
d_f = zeros(fixed_number, 1);

% free decision variable
free_number = r_order*(k_segment-1); % intermedia v,a
d_p = zeros(free_number, 1);

C_t_col = fixed_number + free_number;
C_t_row = k_segment * 2 * (1+r_order); % total number of derivative constraints each segment front and end (p,v,a)

C_t = zeros(C_t_row, C_t_col); % C_t is C^T

c0 = eye(3); % for final and initial p,v,a selection mapping
c1 = zeros(4,1);
c1(1) = 1;
c1(4) = 1; % for intermediate waypoints

c2 = zeros(5,2);
c2(1,1) = 1;
c2(2,2) = 1;
c2(4,1) = 1;
c2(5,2) = 1;
c2_alternative = eye(2);

% for k-segment trajectory
% C contains:
% -> 2 * c0
% -> (k-1) * c1
% -> (k-1) * c2
row_count = 1;
col_count = 1;
C_t(1:3, 1:3) = c0; % fixed position
row_count = row_count + 3;
col_count = col_count + 3;

for i=1:(k_segment - 1)
    C_t(row_count:row_count+3,col_count) = c1;
    row_count = row_count + 6;
    col_count = col_count+1;
end

C_t(row_count:row_count+2, col_count:col_count+2) = c0;

col_count = col_count + 3;
row_count = 5;

for i=1:(k_segment - 1)
    C_t(row_count:row_count+4,col_count:col_count+1) = c2;
    row_count = row_count + 6;
    col_count = col_count+2;
end

C = C_t';

time_vector = allocate_time(waypts, T);
Q_hessian = calQ(time_vector);

M = A;

R = C * inv(M)' * Q_hessian * inv(M) * C_t;

%R_ff (fixed_number*fixed_numer)
R_ff = zeros(fixed_number);
R_fp = zeros(fixed_number, free_number);
R_pf = zeros(free_number, fixed_number);
R_pp = zeros(free_number, free_number);

R_ff = R(1:fixed_number,1:fixed_number);
R_fp = R(1:fixed_number,fixed_number+1:fixed_number+1+free_number-1);
R_pf = R(fixed_number+1:fixed_number+1+free_number-1,1:fixed_number);
R_pp = R(fixed_number+1:fixed_number+1+free_number-1,fixed_number+1:fixed_number+1+free_number-1);

d_f(1) =  waypts(1,1);
for i=4:4+k_segment-1
    d_f(i) = waypts(1,i-2);
end
d_p_primal = -inv(R_pp) * R_fp' * d_f;

d_primal = [d_f;d_p_primal];

d = C_t * d_primal;

p = inv(M) * d






