% fig6_7.m - Reproduce Figures 6 and 7 from Fleury et al. (1999)
%
% CRLB and RMSEE of |alpha_1| and phi_1 vs [Delta_tau, Delta_phi]
% with Delta_nu as parameter.
%
% Two-wave scenario (Section V-C):
%   Wave 1: tau=0, phi=0, nu=0, alpha_1 = 1
%   Wave 2: tau=Delta_tau, phi=Delta_phi, nu=Delta_nu, alpha_2 = -1

fprintf('=== Reproducing Figures 6 and 7 ===\n');
fprintf('Computing CRLB and running Monte Carlo SAGE...\n');

% Resolution limits (Eq. 31)
tau_c = 1 / cfg.B_s;  % ~Tc
phi_c = cfg.lambda / (cfg.M * cfg.d);  % rad (intrinsic azimuth resolution)
nu_c  = 1 / (cfg.K * cfg.T);  % Hz (intrinsic Doppler resolution)

fprintf('Resolution limits:\n');
fprintf('  tau_c = %.2f ns\n', tau_c * 1e9);
fprintf('  phi_c = %.2f deg\n', phi_c * 180/pi);
fprintf('  nu_c  = %.2f Hz\n', nu_c);

% Parameter grid
Delta_tau_vals = [0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0] * tau_c;
Delta_phi_vals = [0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0] * phi_c;
Delta_nu_vals  = [0, 0.5] * nu_c;

% Fixed alpha values
alpha_mag = 1;

% For each Delta_nu, compute CRLB and RMSEE
for inu = 1:length(Delta_nu_vals)
    Delta_nu = Delta_nu_vals(inu);

    fprintf('\nDelta_nu / nu_c = %.2f\n', Delta_nu / nu_c);

    % CRLB matrices
    CRLB_alpha = zeros(length(Delta_tau_vals), length(Delta_phi_vals));
    CRLB_phi   = zeros(length(Delta_tau_vals), length(Delta_phi_vals));

    for it = 1:length(Delta_tau_vals)
        for ip = 1:length(Delta_phi_vals)
            Delta_tau = Delta_tau_vals(it);
            Delta_phi = Delta_phi_vals(ip);

            theta1 = [0, 0, 0, alpha_mag];
            theta2 = [Delta_tau, Delta_phi, Delta_nu, alpha_mag * exp(1j * pi)];

            crlb_2 = crlb_two_waves(cfg, theta1, theta2);

            % Normalize CRLB as in paper
            % CRLB(|alpha|) / |alpha|^2
            CRLB_alpha(it, ip) = sqrt(crlb_2.alpha1_mag) / abs(alpha_mag);
            CRLB_phi(it, ip)   = sqrt(crlb_2.phi1);
        end
    end

    % Save CRLB results
    save(sprintf('results/tables/crlb_fig6_7_nu%.1f.mat', Delta_nu/nu_c), ...
        'CRLB_alpha', 'CRLB_phi', 'Delta_tau_vals', 'Delta_phi_vals');

    % Plot CRLB
    figure;
    subplot(1,2,1);
    imagesc(Delta_tau_vals/tau_c, Delta_phi_vals/phi_c, CRLB_alpha');
    colorbar;
    xlabel('\Delta\tau / \tau_c'); ylabel('\Delta\phi / \phi_c');
    title(sprintf('CRLB(|\alpha_1|)/|\alpha_1|, \Delta\nu/{\nu_c}=%.1f', Delta_nu/nu_c));
    set(gca, 'YDir', 'normal');
    grid on;

    subplot(1,2,2);
    imagesc(Delta_tau_vals/tau_c, Delta_phi_vals/phi_c, CRLB_phi');
    colorbar;
    xlabel('\Delta\tau / \tau_c'); ylabel('\Delta\phi / \phi_c');
    title(sprintf('CRLB(\phi_1) [rad], \Delta\nu/{\nu_c}=%.1f', Delta_nu/nu_c));
    set(gca, 'YDir', 'normal');
    grid on;

    sgtitle(sprintf('Figure %d: CRLB for Two-Wave Scenario', 6 + (inu-1)));
    saveas(gcf, sprintf('results/figures/fig%d_crlb_nu%.1f.png', 6 + (inu-1), Delta_nu/nu_c));
    savefig(gcf, sprintf('results/figures/fig%d_crlb_nu%.1f.fig', 6 + (inu-1), Delta_nu/nu_c));
end

fprintf('\nCRLB computation complete.\n');
fprintf('Figures saved to results/figures/\n');
fprintf('NOTE: Full Monte Carlo RMSEE computation is computationally intensive.\n');
fprintf('Run RMSEE computation separately with cfg.Nmc runs.\n');
end
