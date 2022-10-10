% function to get the mapping between polynomial coefficients (p) and segment
% end point derivative constriants (d)
% A * p = d

function A=getMappingA(n_order, k_segment, t_alloc, dim)
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
end
% 
% A_sub(1,:) = poly_evaluate(0,t_alloc(1),n_order);
% A_sub(2,:) = poly_evaluate(1,t_alloc(1),n_order);
% A_sub(3,:) = poly_evaluate(2,t_alloc(1),n_order);
% A_sub(4,:) = poly_evaluate(0,t_alloc(2),n_order);
% A_sub(5,:) = poly_evaluate(1,t_alloc(2),n_order);
% A_sub(6,:) = poly_evaluate(2,t_alloc(2),n_order);