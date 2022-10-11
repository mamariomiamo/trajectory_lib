%% time allocation function
function t_alloc = allocate_time(waypts,T)

t_alloc = zeros(length(waypts)-1,1);

total_dist = 0; % total distance for all segments

for i=1:(length(waypts)-1)
total_dist = total_dist + norm(waypts(:,i+1)-waypts(:,i));
end

for i=1:length(waypts)-1
    t_alloc(i) = norm(waypts(:,i+1)-waypts(:,i))/total_dist*T;
end
t_alloc = t_alloc';
end