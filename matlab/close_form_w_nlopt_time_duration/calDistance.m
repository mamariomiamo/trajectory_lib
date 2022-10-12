function total_dist = calDistance(waypts, k_segment)
total_dist = 0;
for i=2:(k_segment+1)
    total_dist = total_dist + norm(waypts(:,i) - waypts(:,i-1));
end