% Code for UNSUPERVISED PARTICLE SORTING FOR CRYO-EM USING PROBABILISTIC PCA
% Works on data set of arbitrary dimensions.
% Make sure to have ASPIRE in your path
%
%  Input:   projs     -- size (L,L,N)
%  Output:  est_model -- struct with model parameters
%
% Versions:
% 0.1  |  Gili Weiss-Dicker, December  2020
% 0.2  |  Gili Weiss-Dicker, September 2020 | Add Gram-schmidt method.
% 0.25 |  Gili Weiss-Dicker, March     2021 | Add online sorting.


% Simulation configurations:
% --------------------------
Y = use_steerable_features(projs);                    % Fourier-Bessel coefficents
% problem setting
config.problem_setting.L                = size(Y,1);% High dimension of original data. 
config.problem_setting.K                = 60;       % Low dimension to project data to
config.problem_setting.M                = 1;        % Number of subspaces (Gaussians)
config.problem_setting.N                = size(Y,2);% Number of samples 

% algorithm configurations:
config.algo_setting.maxiter             = 30;       % Maximum EM iterations
config.algo_setting.RandomSeedEMInit    = 0;        % 0 = random seed; otherwise input is seed for EM initialization.
config.algo_setting.GramSchmidt         = 0;        % warning: breaks promise to convergence to local maximum likelihood 
config.algo_setting.Sorting             = 1;
config.algo_setting.SortOut             = 0.05;     % alpha, percentage of images to sort out every sorting step.
config.algo_setting.SortEveryIter       = 6;        % P, number of EM iterations between sorting step.

% running setting
config.run_setting.Verbose              = 1;
config.run_setting.assert               = false;  % must be set to false to run.
  

% EM algorithm
% --------------
[config, est_model] = em_union_of_subspaces(Y, config); 
log_message('DONE!')