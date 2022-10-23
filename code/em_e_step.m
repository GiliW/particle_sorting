function [config, est_model] = em_e_step(Y, model, config)
% EM expectation step
% Output:
%   est_model:   struct. current updated estimated model parameters.

[L, N]              = size(Y);
K                   = size(model.C,2);
M                   = size(model.C,3);
C                   = model.C;
sigma2_gal          = model.sigma2_gal(end);
mu_L                = model.mu_L;
est_model           = model;

E_x_given_yis       = zeros(K, N, M);
E_xxT_given_yis     = zeros(K,K,N,M);
h_i_m               = zeros(N,M);



% Calculate h_i_m
% ---------------
if M > 1
    Sigma    = zeros(L,L,M);
    invSigma = zeros(L,L,M);
    detInvSigma = zeros(1,M);
    for m = 1:M 
        Sigma(:,:,m)         = C(:,:,m)*C(:,:,m).' + sigma2_gal * eye(L);
        invSigma(:,:,m)      = pinv(Sigma(:,:,m)); 
        detInvSigma(m)       = det(invSigma(:,:,m));
    end
    for m = 1:M 
        for ii = 1:N 
            h_i_m(ii,m) = est_model.w_m(m) * detInvSigma(m) * ... 
                exp( -0.5 * ( (Y(:, ii).' - mu_L(:,m).') * ...
                invSigma(:,:,m) * (Y(:, ii).' - mu_L(:,m).').' ) ) + eps;
        end
    end
    h_i_m = h_i_m./ sum(h_i_m,2);
else
   h_i_m = ones(N,M); 
   invSigma = ones(L,L,M);
end

% Estimate E[x|y_i]
% -----------------
B_m = zeros(K,L,M);
for m = 1:M 
    B_m(:,:,m) = C(:,:,m).' * invSigma(:,:,m); 
    if ( sum(find(B_m(:,:,m) == Inf)) || sum(find( isnan(B_m(:,:,m))) ))
        error('ERROR: B is singular, close to singular or badly scaled.')
    end

    E_x_given_yis(:,:,m) = B_m(:,:,m) * (Y - mu_L(:,m)); 
end

% Estimate E[x x^T|y_i]
% ---------------------
for ii = 1:N 
    for m = 1:M 
        E_xxT_given_yis(:,:,ii,m) = eye(K) - B_m(:,:,m) * C(:,:,m) + B_m(:,:,m) * (Y(:, ii) - mu_L(:,m)) * (Y(:, ii) - mu_L(:,m)).'* B_m(:,:,m).'; 
    end
end

% Log-likelihood
% --------------
llh_inner_yi = zeros(N, M);
for ii = 1:N 
    for m = 1:M 

       inner_yi = Y(:,ii);
        llh_inner_yi(ii, m) = h_i_m(ii,m) *  ( inner_yi' * inner_yi - 2 * inner_yi' * mu_L(:,m) + ...
            + 2 * ( mu_L(:,m)' - inner_yi' ) * C(:,:,m) * E_x_given_yis(:,ii,m) + mu_L(:,m)' * mu_L(:,m) + ...
            trace(C(:,:,m)' * C(:,:,m) * E_xxT_given_yis(:,:,ii,m)) );

        if isnan(llh_inner_yi(ii,m))
            error('NaN values in log-likelihood')
        end
    end
end
llh = - 0.5 *  L * log( sigma2_gal +eps) - (1/(2 * N * sigma2_gal + eps))*sum( sum(llh_inner_yi) );
    
est_model.llh             = [est_model.llh, llh];
est_model.E_x_given_yis   = E_x_given_yis;
est_model.E_xxT_given_yis = E_xxT_given_yis; 
est_model.h_i_m           = h_i_m;

end
