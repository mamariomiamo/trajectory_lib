function [R_ff, R_fp, R_pf, R_pp]=getSubR(R, fixed_number, free_number)
R_ff = zeros(fixed_number);
R_fp = zeros(fixed_number, free_number);
R_pf = zeros(free_number, fixed_number);
R_pp = zeros(free_number, free_number);

R_ff = R(1:fixed_number,1:fixed_number);
R_fp = R(1:fixed_number,fixed_number+1:fixed_number+1+free_number-1);
R_pf = R(fixed_number+1:fixed_number+1+free_number-1,1:fixed_number);
R_pp = R(fixed_number+1:fixed_number+1+free_number-1,fixed_number+1:fixed_number+1+free_number-1);
end