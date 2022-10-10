% min sqrt(x2)
% s.t. x2≥0 , x2≥(a1 x1+b1)^3, and x2≥(a2 x1+b2)^3
% a1=2, b1=0, a2=-1, b2=1.

clc, clear;

% opt.algorithm = NLOPT_LD_MMA;
opt.algorithm = NLOPT_LN_COBYLA;

opt.verbose = 1;

opt.lower_bounds = [-inf, 0];

opt.min_objective = @myfunc;

opt.fc = { (@(x) myconstraint(x,2,0)), (@(x) myconstraint(x,-1,1)) };

opt.fc_tol = [1e-8, 1e-8];

opt.xtol_rel = 1e-4;

[xopt, fmin, retcode] = nlopt_optimize(opt, [100 500.678])

[xopt, fmin, retcode] = nlopt_optimize(opt, [0.3 0.3])