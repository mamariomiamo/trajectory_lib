function total_dist = calDistance(waypts, k_segment)
total_dist = 0;
for i=1:k_segment
    total_dist = total_dist + norm(waypts(:,i+1) - waypts(:,i));
end