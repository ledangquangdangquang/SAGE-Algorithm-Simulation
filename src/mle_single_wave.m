function [theta_hat, z_max] = mle_single_wave(x, cfg, search_grid, theta_init)
% [theta_hat, z_max] = mle_single_wave(x, cfg, search_grid, theta_init)
% Maximum Likelihood estimation of a single wave's parameters.
%
% Implements the coordinate-wise updating procedure Eq. (14) with
% parabolic interpolation for sub-grid accuracy (Section IV-A).
% theta_init is the current estimate used as the 1-D search starting point.

M = cfg.M;
K = cfg.K;
s_power = cfg.s_power;

tau_grid = search_grid.tau;
phi_grid = search_grid.phi;
nu_grid  = search_grid.nu;

% Start from current estimate if provided, else first grid point
if nargin >= 4
    tau_hat = theta_init(1);
    phi_hat = theta_init(2);
    nu_hat  = theta_init(3);
else
    tau_hat = tau_grid(1);
    phi_hat = phi_grid(1);
    nu_hat  = nu_grid(1);
end

% ----- 1-D search over tau (with phi_hat, nu_hat from current estimate) -----
z_sq_tau = zeros(size(tau_grid));
for it = 1:length(tau_grid)
    z_val = correlator_output(x, tau_grid(it), phi_hat, nu_hat, cfg);
    z_sq_tau(it) = abs(z_val)^2;
end
[tau_hat, ~] = parabolic_interp(z_sq_tau, tau_grid);

% ----- 1-D search over phi (with refined tau_hat, nu_hat) -----
z_sq_phi = zeros(size(phi_grid));
for ip = 1:length(phi_grid)
    z_val = correlator_output(x, tau_hat, phi_grid(ip), nu_hat, cfg);
    z_sq_phi(ip) = abs(z_val)^2;
end
[phi_hat, ~] = parabolic_interp(z_sq_phi, phi_grid);

% ----- 1-D search over nu (with refined tau_hat, phi_hat) -----
z_sq_nu = zeros(size(nu_grid));
for in = 1:length(nu_grid)
    z_val = correlator_output(x, tau_hat, phi_hat, nu_grid(in), cfg);
    z_sq_nu(in) = abs(z_val)^2;
end
[nu_hat, ~] = parabolic_interp(z_sq_nu, nu_grid);

% Compute amplitude from exact correlation at refined parameters
z_best = correlator_output(x, tau_hat, phi_hat, nu_hat, cfg);
alpha_hat = z_best / (s_power * M * K);

theta_hat = [tau_hat, phi_hat, nu_hat, alpha_hat];
z_max = abs(z_best)^2;
end
