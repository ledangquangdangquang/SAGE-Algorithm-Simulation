function [theta_hist, logl_hist] = sage_algorithm(y, cfg, theta_init, search_grid)
% [theta_hist, logl_hist] = sage_algorithm(y, cfg, theta_init, search_grid)
% Main SAGE algorithm loop.
%
% Implements the SFG in Fig. 5 of Fleury et al. (1999).
% Coordinate-wise updating procedure Eq. (14).
%
% Inputs:
%   y           - M x Nsamp x K received signal
%   cfg         - configuration struct
%   theta_init - L x 4 initial parameter estimates
%   search_grid - struct with fields: tau, phi, nu (grid vectors)
%
% Outputs:
%   theta_hist - (n_iter+1) x L x 4 parameter history
%   logl_hist  - (n_iter+1) x 1 log-likelihood history

L = size(theta_init, 1);
n_iter = cfg.n_iter;

% Preallocate history
theta_hist = zeros(n_iter + 1, L, 4);
theta_hist(1, :, :) = theta_init;

% Log-likelihood history
logl_hist = zeros(n_iter + 1, 1);
logl_hist(1) = log_likelihood(y, theta_init, cfg);

theta = theta_init;

for iter = 1:n_iter
    for l = 1:L
        % E-step: compute complete data estimate for wave l
        X_hat_l = expectation_step(y, theta, l, cfg);

        % M-step: estimate wave l parameters (start from current estimate)
        theta_l_new = mle_single_wave(X_hat_l, cfg, search_grid, theta(l, :));

        % Update
        theta(l, :) = theta_l_new;
    end

    % Store history
    theta_hist(iter + 1, :, :) = theta;

    % Compute log-likelihood
    logl_hist(iter + 1) = log_likelihood(y, theta, cfg);
end
end
