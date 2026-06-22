function [y, x_true] = signal_model(cfg, theta_true)
% [y, x_true] = signal_model(cfg, theta_true)
% Generate received signal according to Eq. (2).
%
% y(t) = sum_{l=1}^L x_l(t; theta_l) + sqrt(N0/2) * w(t)
%
% Inputs:
%   cfg        - configuration struct (from config.m)
%   theta_true - L x 4 matrix [tau, phi, nu, alpha] for each wave
%
% Outputs:
%   y     - M x N_samp x K array of received signals
%   x_true - M x N_samp x K x L array of individual wave contributions

M = cfg.M;
K = cfg.K;
P = cfg.P;
Tc = cfg.Tc;
Ts = cfg.Ts;
Ns = cfg.Ns;
d = cfg.d;
lambda = cfg.lambda;
s = cfg.s;
s_power = cfg.s_power;

% Allocate
y = zeros(M, Ns, K);
x_true = zeros(M, Ns, K, size(theta_true, 1));

% Time vector for one burst
t = (0:Ns-1).' * Ts;

for l = 1:size(theta_true, 1)
    tau_l = theta_true(l, 1);
    phi_l = theta_true(l, 2);
    nu_l  = theta_true(l, 3);
    alpha_l = theta_true(l, 4);

    % Steering vector Eq. (1)
    c_phi = steering_vector(phi_l, M, d, lambda);

    for k = 1:K
        % Time base for this snapshot
        t_center = (k - 1) * P * Tc;
        t_k = t + t_center;

        % Doppler phase
        doppler = exp(1j * 2 * pi * nu_l * t_k);

        % Delayed sounding sequence
        tau_samp = round(tau_l / Ts);
        s_delayed = circshift(s, tau_samp);

        % Wave contribution Eq. (1)
        x_k = alpha_l * (c_phi * s_delayed.') .* (doppler.');
        x_true(:, :, k, l) = x_k;
    end
end

% Sum all waves
x_sum = sum(x_true, 4);

% Compute noise variance from SNR definition
% SNR per wave per antenna = |alpha_l|^2 * ||s||^2 / (N0/2)
% N0 = 2 * |alpha_l|^2 * ||s||^2 / SNR_lin
N0 = 2 * abs(theta_true(1, 4))^2 * s_power / cfg.SNR_lin;

% Generate noise (use rng_seed from cfg, or incrementing counter if provided)
if isfield(cfg, 'mc_run')
    noise_seed = cfg.rng_seed + cfg.mc_run;
else
    noise_seed = cfg.rng_seed + 1000;
end
rng(noise_seed);
noise = sqrt(N0 / 2) * (randn(M, Ns, K) + 1j * randn(M, Ns, K));

% Received signal Eq. (2)
y = x_sum + noise;
end
