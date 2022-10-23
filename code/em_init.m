function est_model = em_init(Y, K, M, config)
% Random initialization: which Gaussian does each sample belong to.
% Input: 
%   Y:       L x n data matrix
%   K:       1 x 1 dimension to project data to
%   M:       number of Gaussians
% Output:
%   est_model:  random initialization for model parameters.


if config.algo_setting.RandomSeedEMInit == 0 % not reproducible
    rng('shuffle', 'twister'); 
else % seed
    rng(config.algo_setting.RandomSeedEMInit);   
end
L                                = size(Y,1);
est_model.K                      = K;
est_model.Y                      = Y;
est_model.inliers_ind            = 1:size(Y,2);
est_model.subspace_iter          = 1;

% model parameters
% -----------------
% est_model = {};
C = zeros(size(Y,1), K, M);
W0 = sparse(eye(L));   % hyperparameter of inverse Wishart prior of covariances
v0 = L+1;              % hyperparameter of inverse Wishart prior of covariances
for m = 1:M
    PSD_mat = iwishrnd(W0,v0);
    C(:,:,m) = normc(PSD_mat(1:L,1:K)); % normalize C's columns
end
est_model.C          = C;
% est_model.C                      = randn(L,K,M); % subspaces matrices
est_model.mu_L                   = randn(L,M);   % subspaces means
est_model.E_x_given_yis          = zeros(1,K,M);
est_model.E_xxT_given_yis        = zeros(K,K,M); 
for m = 1:M
    est_model.E_xxT_given_yis(:,:,m) = eye(K,K); 
end
est_model.w_m                    = ones(1,M)/M; % mixing weights

% model noise
% -----------
est_model.sigma2_gal = abs(randn(1));      % noise variance

% log-likelihood 
% --------------
est_model.llh        = -Inf;   

if config.run_setting.Verbose
    log_message('EM: model initialization done.');
end

end