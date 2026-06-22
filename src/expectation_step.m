function X_hat_l = expectation_step(y, theta_hat_all, l, cfg)
% X_hat_l = expectation_step(y, theta_hat_all, l, cfg)
% Compute the conditional expectation (E-step) for wave l.
%
% Eq. (13) in Fleury et al. (1999)
% X_hat_l = x_l(theta_hat_l) + beta_l * (y - sum_{i=1}^L x_i(theta_hat_i))
%
% Inputs:
%   y              - M x Nsamp x K observed signal
%   theta_hat_all  - L x 4 matrix of current parameter estimates
%   l              - index of wave to re-estimate
%   cfg            - configuration struct
%
% Output:
%   X_hat_l - M x Nsamp x K complete data estimate for wave l

M = cfg.M;
K = cfg.K;
P = cfg.P;
Tc = cfg.Tc;
Ts = cfg.Ts;
beta_l = cfg.beta;
Ns = cfg.Ns;
d = cfg.d;
lambda = cfg.lambda;
s = cfg.s;

% Time vector
t = (0:Ns-1).' * Ts;

% Compute contribution of each wave
L = size(theta_hat_all, 1);
x_hat_all = zeros(M, Ns, K, L);

for i = 1:L
    tau_i = theta_hat_all(i, 1);
    phi_i = theta_hat_all(i, 2);
    nu_i  = theta_hat_all(i, 3);
    alpha_i = theta_hat_all(i, 4);

    c_phi_i = steering_vector(phi_i, M, d, lambda);

    for k = 1:K
        t_center = (k - 1) * P * Tc;
        t_k = t + t_center;
        doppler = exp(1j * 2 * pi * nu_i * t_k);
        tau_samp = round(tau_i / Ts);
        s_delayed = circshift(s, tau_samp);
        x_hat_all(:, :, k, i) = alpha_i * (c_phi_i * s_delayed.') .* (doppler.');
    end
end

% Sum of all current estimates
x_sum = sum(x_hat_all, 4);

% Eq. (13): X_hat_l = x_l + beta_l * (y - sum_i x_i)
X_hat_l = x_hat_all(:, :, :, l) + beta_l * (y - x_sum);
end
