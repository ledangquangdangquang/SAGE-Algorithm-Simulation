function crlb = crlb_single_wave(cfg, theta)
% crlb = crlb_single_wave(cfg, theta)
% Compute CRLB for a single wave's parameters.
%
% Eqs. (27)-(30) in Fleury et al. (1999):
%   CRLB(tau)  = 1/(8*pi^2*rho*B_g^2)          (27)
%   CRLB(phi)  = 1/(2*rho*M*(pi*d/lambda)^2*(M^2-1)/3)  (28)
%   CRLB(nu)   = 3/(2*pi^2*rho*K*Ts^2*(K^2-1))        (29)
%   CRLB(|alpha|) = |alpha|^2/(2*rho)                    (30)
%
% where rho = |alpha|^2 * ||s||^2 * M * K / N0 (SNR at correlator output)
%
% Inputs:
%   cfg   - configuration struct
%   theta - [tau, phi, nu, alpha] for the wave
%
% Output:
%   crlb - struct with fields: tau, phi, nu, alpha_mag

M = cfg.M;
K = cfg.K;
P = cfg.P;
Tc = cfg.Tc;
Ts = cfg.Ts;
d = cfg.d;
lambda = cfg.lambda;
s = cfg.s;
s_power = cfg.s_power;

% Extract parameters
tau   = theta(1);
phi   = theta(2);
nu    = theta(3);
alpha = theta(4);

% SNR at correlator output (rho in paper)
% rho = |alpha|^2 * ||s||^2 * M * K / N0
% N0 = 2 * |alpha|^2 * s_power / SNR_lin
N0 = 2 * abs(alpha)^2 * s_power / cfg.SNR_lin;
rho = abs(alpha)^2 * s_power * M * K / N0;

% Gabor bandwidth of s(t) (computed numerically)
% B_g^2 = integral f^2 |S(f)|^2 df / integral |S(f)|^2 df
Nfft = 2^nextpow2(length(s));
S = fft(s, Nfft);
f = ((0:Nfft-1)/Nfft * (1/Ts)).';
f_centered = f - 1/(2*Ts);
S_mag_sq = abs(S).^2;
B_g_sq = sum(f_centered.^2 .* S_mag_sq) / sum(S_mag_sq);
B_g = sqrt(B_g_sq);

% CRLB(tau) - Eq. (27)
crlb.tau = 1 / (8 * pi^2 * rho * B_g_sq);

% CRLB(phi) - Eq. (28)
% For ULA with reference at center:
% F_pp = 2*rho*M*(pi*d/lambda)^2 * (M^2-1)/3
% CRLB(phi) = 1/F_pp
F_pp = 2 * rho * M * (pi * d / lambda * cos(phi))^2 * (M^2 - 1) / 12;
crlb.phi = 1 / F_pp;

% CRLB(nu) - Eq. (29)
crlb.nu = 3 / (2 * pi^2 * rho * K * Ts^2 * (K^2 - 1));

% CRLB(|alpha|) - Eq. (30)
crlb.alpha_mag = abs(alpha)^2 / (2 * rho);

% RMSEE (root)
crlb.rmsee_tau = sqrt(crlb.tau);
crlb.rmsee_phi = sqrt(crlb.phi);
crlb.rmsee_nu  = sqrt(crlb.nu);
crlb.rmsee_alpha = sqrt(crlb.alpha_mag);
end
