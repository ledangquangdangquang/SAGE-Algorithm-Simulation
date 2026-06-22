% fig10_11.m - Reproduce Figures 10 and 11 from Fleury et al. (1999)
%
% RMSEE vs number of iteration cycles for different [Delta_tau, Delta_phi]
% and Delta_nu values.
%
% Fig 10: Delta_nu = 0
% Fig 11: Delta_nu = 0.5 * nu_c

fprintf('=== Reproducing Figures 10 and 11 ===\n');

% Resolution limits
tau_c = 1 / cfg.B_s;
phi_c = cfg.lambda / (cfg.M * cfg.d);
nu_c  = 1 / (cfg.K * cfg.T);

fprintf('Resolution limits:\n');
fprintf('  tau_c = %.2f ns\n', tau_c * 1e9);
fprintf('  phi_c = %.2f deg\n', phi_c * 180/pi);
fprintf('  nu_c  = %.2f Hz\n', nu_c);

% Parameter pairs [Delta_tau/tau_c, Delta_phi/phi_c] to test
param_pairs = [
    0.5, 0;    % only delay separation
    1.0, 0;
    0,   0.5;  % only azimuth separation
    0,   1.0;
    0.5, 0.5;  % both
    1.0, 1.0;
];

% Doppler separations for Fig 10 (nu=0) and Fig 11 (nu=0.5*nu_c)
nu_seps = [0, 0.5 * nu_c];

% Search grid (limited range)
grd = compute_search_grid(cfg);
grd.tau = 0:cfg.dtau:(3 * tau_c);
grd.phi = (-5*phi_c):cfg.dphi:(5*phi_c);
grd.nu = (-nu_c):cfg.dnu:(nu_c);

cfg.n_iter = 30;  % Enough for convergence plot

for inu = 1:length(nu_seps)
    Delta_nu = nu_seps(inu);
    fig_num = 9 + inu;  % Fig 10 or 11

    fprintf('\n=== Delta_nu / nu_c = %.2f (Figure %d) ===\n', ...
        Delta_nu / nu_c, fig_num);

    % For each parameter pair, run Monte Carlo and compute RMSEE history
    rmsee_alpha_hist = zeros(size(param_pairs, 1), cfg.n_iter + 1);
    rmsee_phi_hist   = zeros(size(param_pairs, 1), cfg.n_iter + 1);

    for ip = 1:size(param_pairs, 1)
        dt_norm = param_pairs(ip, 1);
        dp_norm = param_pairs(ip, 2);

        Delta_tau = dt_norm * tau_c;
        Delta_phi = dp_norm * phi_c;

        fprintf('  [Delta_tau/tau_c=%.2f, Delta_phi/phi_c=%.2f]...', dt_norm, dp_norm);

        % True parameters
        theta1 = [0, 0, 0, 1];
        theta2 = [Delta_tau, Delta_phi, Delta_nu, exp(1j * pi)];

        % Accumulators for RMSEE
        alpha_err_sq = zeros(1, cfg.n_iter + 1);
        phi_err_sq   = zeros(1, cfg.n_iter + 1);
        mc_count = 0;

        for mc = 1:cfg.Nmc
            cfg.mc_run = mc;
            % Generate received signal
            [y, ~] = signal_model(cfg, [theta1; theta2]);

            % Initialize
            theta_init = init_successive_cancellation(y, cfg, grd);

            % Run SAGE
            [theta_hist, ~] = sage_algorithm(y, cfg, theta_init, grd);

            % Compute squared errors for wave 1
            for it = 1:cfg.n_iter + 1
                theta_est = squeeze(theta_hist(it, 1, :)).';

                % Delay and Doppler should converge to theta1
                % But the assignment of waves may swap, so check both
                err_alpha1 = abs(theta_est(4)) - abs(theta1(4));
                err_phi1 = angle_diff(theta_est(2), theta1(2));

                alpha_err_sq(it) = alpha_err_sq(it) + err_alpha1^2;
                phi_err_sq(it)   = phi_err_sq(it)   + err_phi1^2;
            end
            mc_count = mc_count + 1;
        end

        % RMSEE
        rmsee_alpha_hist(ip, :) = sqrt(alpha_err_sq / mc_count);
        rmsee_phi_hist(ip, :)   = sqrt(phi_err_sq / mc_count);
        fprintf(' done.\n');
    end

    % Save results
    save(sprintf('results/tables/rmsee_fig%d.mat', fig_num), ...
        'rmsee_alpha_hist', 'rmsee_phi_hist', 'param_pairs');

    % Plot Fig 10 or 11
    figure;
    subplot(1,2,1);
    colors = lines(size(param_pairs, 1));
    for ip = 1:size(param_pairs, 1)
        dt_norm = param_pairs(ip, 1);
        dp_norm = param_pairs(ip, 2);
        semilogy(0:cfg.n_iter, rmsee_alpha_hist(ip, :), ...
            'o-', 'Color', colors(ip,:), 'LineWidth', 1.2, ...
            'DisplayName', sprintf('\Delta\\tau/\\tau_c=%.1f, \Delta\\phi/\\phi_c=%.1f', dt_norm, dp_norm));
        hold on;
    end
    xlabel('Iteration cycle'); ylabel('RMSEE(|\alpha_1|)/|\alpha_1|');
    legend('Location', 'northeast'); grid on;
    title(sprintf('Fig %d(a): Amplitude convergence', fig_num));

    subplot(1,2,2);
    for ip = 1:size(param_pairs, 1)
        dt_norm = param_pairs(ip, 1);
        dp_norm = param_pairs(ip, 2);
        semilogy(0:cfg.n_iter, rmsee_phi_hist(ip, :) * 180/pi, ...
            'o-', 'Color', colors(ip,:), 'LineWidth', 1.2, ...
            'DisplayName', sprintf('\Delta\\tau/\\tau_c=%.1f, \Delta\\phi/\\phi_c=%.1f', dt_norm, dp_norm));
        hold on;
    end
    xlabel('Iteration cycle'); ylabel('RMSEE(\phi_1) [deg]');
    legend('Location', 'northeast'); grid on;
    title(sprintf('Fig %d(b): Azimuth convergence', fig_num));

    sgtitle(sprintf('Figure %d: SAGE Convergence, \\Delta\\nu/\\nu_c = %.2f', fig_num, Delta_nu/nu_c));
    saveas(gcf, sprintf('results/figures/fig%d.png', fig_num));
    savefig(gcf, sprintf('results/figures/fig%d.fig', fig_num));
end

fprintf('\nFigures 10 and 11 saved to results/figures/\n');
end

function d = angle_diff(a, b)
    d = angle(exp(1j*(a-b)));
end
