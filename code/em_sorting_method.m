function [model, config] = em_sorting_method(model, config, percentToSort) 
% Sorting angel method
sort_score = zeros(1, config.problem_setting.N);
for m = 1:config.problem_setting.M
    x_hat  = ( model.Y.' * gram_schmidt(   model.C(:,:,m) ) ).'  ;
    y_hat = x_hat.' * gram_schmidt( model.C(:,:,m)) .';

    projection_error = vecnorm(y_hat.' - model.Y);
    sort_score = sort_score + projection_error./ vecnorm(x_hat);
end
    
N_top_particles_to_take = floor(config.problem_setting.N*(1-percentToSort));
[~,ii_sorted] = sort(sort_score);

est_inliers_ind  = sort(ii_sorted (1: N_top_particles_to_take));
est_outliers_ind = sort(ii_sorted (N_top_particles_to_take+1 : end));

if ( N_top_particles_to_take > 0 )
    model.inliers_ind     = est_inliers_ind;
    model.est_outliers_ind = est_outliers_ind;

    model.Y               = model.Y(:, est_inliers_ind);
    config.problem_setting.N  = size(model.Y,2);
    model.E_x_given_yis   = model.E_x_given_yis(:,est_inliers_ind,:); 
    model.E_xxT_given_yis = model.E_xxT_given_yis(:,:,est_inliers_ind,:);
    model.h_i_m           = model.h_i_m(est_inliers_ind,:);
end

end
