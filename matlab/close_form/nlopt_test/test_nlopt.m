% min sqrt(x2)
% s.t. x2≥0 , x2≥(a1 x1+b1)^3, and x2≥(a2 x1+b2)^3
% constraints in the form myconstraint(x) ≤ 0
% a1=2, b1=0, a2=-1, b2=1.
% in the original example, myconstraint function is used to evaluate
% nonlinear constraints
clc, clear;

% opt.algorithm = NLOPT_LD_MMA;
opt.algorithm = NLOPT_LN_COBYLA;

opt.verbose = 1;

opt.lower_bounds = [-inf, 0];

f = @(x) sqrt(x);
a = 10;
b = 20;
% myfunc = @(x) (sqrt(x(2)) + a - b);

opt.min_objective = @(x) myfunc(x,a,b);
% opt.min_objective = sqrt(x(2));

opt.fc = { (@(y) myconstraint(y,2,0)), (@(y) myconstraint(y,-1,1)) };

opt.fc_tol = [1e-8, 1e-8];

opt.xtol_rel = 1e-4;

[xopt, fmin, retcode] = nlopt_optimize(opt, [100 500.678])

[xopt, fmin, retcode] = nlopt_optimize(opt, [0.3 0.3])