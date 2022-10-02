clc;clear all;
syms t
r_derivative = 4;
n_order = 5;

q_row = sym(zeros((n_order+1), 1)');

for i=sym((r_derivative+1):(n_order+1))
    q_row(i) = factorial(i-1)/factorial(i-1-r_derivative) * t^(i-1-r_derivative);
end

t = 2;

subs(q_row' * q_row)

int(q_row' * q_row)

q_row' * q_row
