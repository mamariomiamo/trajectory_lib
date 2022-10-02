clc;clear all;
syms t t0 t1
r_derivative = 4;
n_order = 5;

q_row = sym(zeros((n_order+1), 1)');
q_hat_row = sym(zeros((n_order+1), 1)');
q_hat = sym(zeros((n_order+1), (n_order+1)));
q_hat_int = sym(zeros((n_order+1), (n_order+1)));



for i=sym((r_derivative+1):(n_order+1))
    q_row(i) = factorial(i-1)/factorial(i-1-r_derivative) * t^(i-1-r_derivative);
end

for i=sym((r_derivative+1):(n_order+1))
    for j=sym((r_derivative+1):(n_order+1))
        q_hat(i,j) = factorial(i-1)/factorial(i-1-r_derivative) * t^(i-1-r_derivative) * factorial(j-1)/factorial(j-1-r_derivative) * t^(j-1-r_derivative);
    end

end

for i=sym((r_derivative+1):(n_order+1))
    for j=sym((r_derivative+1):(n_order+1))
        q_hat_int(i,j) = factorial(i-1)/factorial(i-1-r_derivative) * factorial(j-1)/factorial(j-1-r_derivative)/(j-1-r_derivative + i-1-r_derivative + 1) * t^(j-1-r_derivative + i-1-r_derivative + 1);
    end
end

% q_hat = q_row.' * q_row
q_hat
q_hat_matlab = q_row.' * q_row
q_hat_int_eval_matlab = int(q_hat, t0, t1)
% q_hat_int
% q_hat_int
t = 0;
a = subs(q_hat_int);

t = 1.05079;
b = subs(q_hat_int);
q_hat_int_eval = vpa(b-a)
