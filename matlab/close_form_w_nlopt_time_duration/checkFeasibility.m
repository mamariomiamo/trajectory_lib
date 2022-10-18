function derivative_initial_guess = checkFeasibility(derivative_initial_guess,v_max_abs, a_max_abs,k_segment)
check_iter = k_segment - 1;
for i=1:length(derivative_initial_guess)
    if abs(derivative_initial_guess(i)) > v_max_abs
        derivative_initial_guess(i) = sign(derivative_initial_guess(i)) *v_max_abs;
    end
end
