function     [C_m_final, mu_L_final, config] = em_linear_eq_solver(config, Y,M,N, E_x_given_yis, E_xxT_given_yis, model, C_m_0, mu_L_0)
% Calculate C_M, mu_L iteratively.
% Input: 
%   Y:                      L x N data matrix
%   M:                      number of Gaussians
%   N:                      number of samples
%   E_x_given_yis:          expected x for each y      
%   E_xxT_given_yis:        expected xx^T for each y    
%   model:                  current model parameters (h_i_m, K)
%   C_m_0:                  init value for C_m            
%   mu_L_0:                 init value for mu_L
%
% Output:
%   C_m:                    new estimation for C_m
%   mu_L:                   new estimation for mu_L
% L        = size(Y,1);
K        = size(C_m_0,2);
C_m_final = zeros(size(C_m_0)); mu_L_final = zeros(size(mu_L_0));

 
for m = 1:M 
    sum_E_xxT_given_yis_m       = reshape( sum( model.h_i_m(:,m) .* reshape( E_xxT_given_yis(:,:,:,m) ,[] ,N)')', K,K) ;
    pinv_sum_E_xxT_given_yis_m  = pinv(sum_E_xxT_given_yis_m );
    
    a = 0; b = 0; a_2 = 0; b_2 = 0;
    for i =1:N
        a = a + (model.h_i_m(i,m).' * Y(:,i).').' *  E_x_given_yis(:,i).';
        b = b + model.h_i_m(i,m) *  E_x_given_yis(:,i).';

        a_2 = a_2 + (model.h_i_m(i,m).' * Y(:,i).').' ;
        b_2 = b_2 + model.h_i_m(i,m) *  E_x_given_yis(:,i);
    end
    a = a * pinv_sum_E_xxT_given_yis_m;
    b = b * pinv_sum_E_xxT_given_yis_m;

    
    mu_L_final(:,m)  = (a_2 - a*b_2)/(N - b*b_2);
    C_m_final(:,:,m) = (a - mu_L_final(:,m) * b ) * pinv_sum_E_xxT_given_yis_m; 
end

end