function time = plot_time(time_vector)
time = zeros(1,length(time_vector)+1);
for i=2:length(time)
    time(i) = time(i-1) + time_vector(i-1);
end