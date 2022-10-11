function C=getSelectionC(fixed_number, free_number, k_segment, r_order)
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
% c2_alternative = eye(2);

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
end