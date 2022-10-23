function [config, est_model] = em_m_step(Y, model, config)
% EM maximization step
% Output:
%   est_model: struct. current updated estimated model parameters

[L, N]          = size(Y);
M               = size(model.C,3);
E_x_given_yis   = model.E_x_given_yis;
E_xxT_given_yis = model.E_xxT_given_yis;
est_model       = model;

sigma2_gal_yi_m = zeros(N,M);

% Solve linear equations for C and mu
% -----------------------------------
[C_m, mu_L, config] = em_linear_eq_solver(config,Y,M,N, E_x_given_yis, E_xxT_given_yis, model, model.C, model.mu_L);


% Calculate model noise
% ----------------------
for ii = 1:N
    for m = 1:M 
        sigma2_gal_yi_m(ii,m) = sigma2_gal_yi_m(ii,m) + ...
             model.h_i_m(ii,m) * ...
             ( vecnorm(Y(:,ii) - mu_L(:,m) ) + ...
               trace(  C_m(:,:,m).' *  C_m(:,:,m)*E_xxT_given_yis(:,:,ii,m) ) +... 
               2 * (mu_L(:,m) - Y(:,ii) ).' *  C_m(:,:,m) * E_x_given_yis(:,ii,m) ); 
    end
end
sigma2_gal = (1/L)* sum(sum(sigma2_gal_yi_m)); 

% Calculate w_m
% -------------
est_model.w_m = sum( model.h_i_m, 1)/N;
est_model.w_m = est_model.w_m./ sum(est_model.w_m);

est_model.C          = C_m; 
est_model.mu_L       = mu_L; 
est_model.sigma2_gal = [est_model.sigma2_gal, sigma2_gal];
end