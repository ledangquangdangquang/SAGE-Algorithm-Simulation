function grid = compute_search_grid(cfg)
% grid = compute_search_grid(cfg)
% Compute parameter search grids for SAGE algorithm.
% Based on quantization precision in Section V-D.

P = cfg.P;
Tc = cfg.Tc;
Ts = cfg.Ts;
K = cfg.K;

% Delay grid: [0, P*Tc] with step dtau
grid.tau = 0:cfg.dtau:(P * Tc);

% Azimuth grid: [-pi/2, pi/2] with step dphi
grid.phi = (-pi/2):cfg.dphi:(pi/2);

% Doppler grid: [-nu_max, nu_max] with step dnu
% Maximum Doppler: 1/(2*Ts) from Nyquist
nu_max = 1 / (2 * Ts);
grid.nu = (-nu_max):cfg.dnu:(nu_max);
end
