%% Calculate matrix Q
% Q is the hessian for min_snap of a 5th order polynomial 
function Q_total = getHessian(t_alloc)
Q = zeros(6,6,length(t_alloc));
Q_total = [];
for i=1:length(t_alloc)
    Q(5,5,i) = 24^2*(t_alloc(i)-0);
    Q(5,6,i) = 0.5 * 24 * 120 * (t_alloc(i)^2 - 0^2);
    Q(6,5,i) = Q(5,6,i);
    Q(6,6,i) = 1/3*120^2*(t_alloc(i)^3 - 0^3);
    Q_total = blkdiag(Q_total,Q(:,:,i));
end
end