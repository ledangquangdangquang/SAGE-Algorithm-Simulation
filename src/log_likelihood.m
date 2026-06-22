function logl = log_likelihood(y, theta, cfg)
% logl = log_likelihood(y, theta, cfg)
% Compute log-likelihood function Lambda(theta; y) from Eq. (4).
%
% Lambda(theta; y) = - (1/N0) * || y - sum_l x_l(theta_l) ||^2 + constant
%
% Inputs:
%   y     - M x Nsamp x K received signal
%   theta - L x 4 parameter matrix
%   cfg   - configuration struct
%
% Output:
%   logl  - scalar log-likelihood value

M = cfg.M;
K = cfg.K;
P = cfg.P;
Tc = cfg.Tc;
Ts = cfg.Ts;
Ns = cfg.Ns;
d = cfg.d;
lambda = cfg.lambda;
s = cfg.s;

% Reconstruct signal from current estimates
x_recon = zeros(M, Ns, K);
t = (0:Ns-1).' * Ts;

for i = 1:size(theta, 1)
    tau_i = theta(i, 1);
    phi_i = theta(i, 2);
    nu_i  = theta(i, 3);
    alpha_i = theta(i, 4);

    c_phi_i = steering_vector(phi_i, M, d, lambda);

    for k = 1:K
        t_center = (k - 1) * P * Tc;
        t_k = t + t_center;
        doppler = exp(1j * 2 * pi * nu_i * t_k);
        tau_samp = round(tau_i / Ts);
        s_delayed = circshift(s, tau_samp);
        x_recon(:, :, k) = x_recon(:, :, k) + ...
            alpha_i * (c_phi_i * s_delayed.') .* (doppler.');
    end
end

% Compute noise variance (same as in signal_model)
N0 = 2 * abs(theta(1, 4))^2 * cfg.s_power / cfg.SNR_lin;

% Residual energy
residual = y - x_recon;
residual_norm = sum(abs(residual(:)).^2);

% Log-likelihood (up to additive constant)
logl = -residual_norm / N0;
end
