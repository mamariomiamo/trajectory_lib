function d_primal = getDerivativePrimal_3axis(d_f, R_pp, R_fp, k_segment, waypts)

d_f(1,:) =  waypts(:,1);
for i=4:4+k_segment-1
    d_f(i,:) = waypts(:,i-2);
end
d_p_primal = -inv(R_pp) * R_fp' * d_f;

d_primal = [d_f;d_p_primal];
end