function options = SBLSet()
% Set default options for SBL_v3p12.m

% convergence threshold
options.convergence.error   = 10^(-4);

options.convergence.delay   = 200;

% maximum number of iterations
% if flag == 1, code issues a warning if not convereged to
% error specified by report.convergence
options.convergence.maxiter = 1000;

% solution only accepted after this iteration
options.convergence.min_iteration = 15; 

% status report every xx iterations
options.status_report = 150;

% noise power initialization guess
options.noisepower.guess = 0.1;

% fixed point acceleration [0.5 2]
options.fixedpoint = 2;

% flag for issuing warnings
options.flag = 0;

options.tic = 1;
