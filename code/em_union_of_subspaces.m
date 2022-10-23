% This function optimizes PPCA model parameters by expectation maximization (EM) algorithm.
% Works on data set of arbitrary dimensions.
%
% Input: 
%   Y:                      L x N data matrix
%   config:                 algorithm and model configurations
%   model:                  if present - uses model parameters for initialization ("easy" mode).
%
% Output:
%   est_model:              trained model structure
%
% Versions:
% 0.1 |  Gili Weiss-Dicker, July 2020

function [config, est_model] = em_union_of_subspaces(Y, config)
    
    if config.run_setting.Verbose
        log_message('EM: running ... \n');
    end
    tic;
     
    % Check input
    % -----------
    if isfield(config.algo_setting, 'maxiter')
        max_iter = config.algo_setting.maxiter;
    else
        max_iter = 500; % default value
    end
    if ~isfield(config.algo_setting, 'SubsetsForEM')
        config.algo_setting.SubsetsForEM = 1; % use full dataset for each EM iter.
    end
    
    K = config.problem_setting.K; % low dimension to project data to
    M = config.problem_setting.M; % number of Gaussians
    
    % Init
    % ----
    tol       = 1e-3;
    est_model = em_init(Y, K, M, config); 
    
    
    % Run
    % ---
    for iter = 1:max_iter
        config.run_setting.iter = iter;
        
        % E step
        % ------
        [config, est_model] = em_e_step(est_model.Y, est_model, config); 
        if config.run_setting.Verbose
            log_message('EM: expectation  step done.');
        end
        if config.run_setting.assert == true
            log_message('ASSERTION of EM algortihm! Run into NaN.1');
            return;
        end
        
        % M step
        % ------
        [config, est_model] = em_m_step(est_model.Y, est_model, config); 
        if config.run_setting.Verbose
            log_message('EM: maximization step done.');
        end

        % Gram Schmidt
        % ------------
        if config.algo_setting.GramSchmidt
            for m = 1:M
                est_model.C(:,:,m) = GramSchmidt(est_model.C(:,:,m)); 
            end
        end

        % Sorting
        % -------
        if config.algo_setting.Sorting && config.run_setting.iter > 1 && mod(config.run_setting.iter,config.algo_setting.SortEveryIter)==0
            log_message('EM: Sorting!')
            [est_model, config] = em_sorting_method(est_model, config, config.algo_setting.SortOut); 
        end
    
        if iter>1 && abs(est_model.llh(iter)-est_model.llh(iter-1)) < tol*abs(est_model.llh(iter)) 
            log_message('EM: Sorting done!')
            break; 
        end
    end
    est_model.timeTook = toc;
end