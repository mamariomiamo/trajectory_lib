clc; close all;

H = [1 -1; -1 2]; 
f = [-2; -6];
A = [1 1; -1 2; 2 1];
b = [2; 2; 3];
% [x, fval] = quadprog(q,p,A,b,Aeq,beq,[1;0],[inf;0]);
[x,fval,exitflag,output,lambda] =quadprog(H,f,A,b)
