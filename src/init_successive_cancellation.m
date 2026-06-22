function theta_init = init_successive_cancellation(y, cfg, search_grid)
% theta_init = init_successive_cancellation(y, cfg, search_grid)
% Initialize SAGE algorithm using successive interference cancellation.
%
% Section IV-B, Method 1 in Fleury et al. (1999).
%
% During the initialization cycle, noncoherent estimation is used:
% no phase information, so we use Eqs. (16)-(17).

L = cfg.L;
M = cfg.M;
P = cfg.P;
Tc = cfg.Tc;
Ts = cfg.Ts;
Ns = cfg.Ns;
K = cfg.K;
d = cfg.d;
lambda = cfg.lambda;
s = cfg.s;
s_power = cfg.s_power;

theta_init = zeros(L, 4);

% Residual signal starts as the observed signal
y_res = y;

for l = 1:L
    % Noncoherent estimation: estimate |z|^2 over tau and phi only
    % (Doppler assumed zero for initialization, Eq. 16)
    tau_grid = search_grid.tau;
    phi_grid = search_grid.phi;

    z_sq_max = -inf;
    tau_hat = 0;
    phi_hat = 0;
    z_at_max = 0;

    for it = 1:length(tau_grid)
        for ip = 1:length(phi_grid)
            z_val = correlator_output(y_res, tau_grid(it), phi_grid(ip), 0, cfg);
            z_sq = abs(z_val)^2;
            if z_sq > z_sq_max
                z_sq_max = z_sq;
                z_at_max = z_val;
                tau_hat = tau_grid(it);
                phi_hat = phi_grid(ip);
            end
        end
    end

    % Estimate amplitude (coherent after location determined)
    alpha_hat = z_at_max / (s_power * M * K);

    theta_init(l, :) = [tau_hat, phi_hat, 0, alpha_hat];

    % Subtract estimated wave from residual (interference cancellation)
    c_phi = steering_vector(phi_hat, M, d, lambda);
    t = (0:Ns-1).' * Ts;

    for k = 1:K
        t_center = (k - 1) * P * Tc;
        t_k = t + t_center;
        tau_samp = round(tau_hat / Ts);
        s_delayed = circshift(s, tau_samp);
        x_l = alpha_hat * (c_phi * s_delayed.') .* (exp(1j * 2 * pi * 0 * t_k).');
        y_res(:, :, k) = y_res(:, :, k) - x_l;
    end
end
end
