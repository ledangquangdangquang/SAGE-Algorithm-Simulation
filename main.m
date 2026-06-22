% main.m - Entry point for SAGE algorithm simulation
% Reproduces results from Fleury et al. (1999)
%
% Usage:
%   main              - run default two-wave scenario demo
%   main('fig6')      - reproduce Fig 6 (CRLB/RMSEE vs [dt, dp])
%   main('fig7')      - reproduce Fig 7
%   main('fig10')     - reproduce Fig 10 (convergence)
%   main('fig11')     - reproduce Fig 11 (convergence)

function main(fig_name)
if nargin < 1
    fig_name = 'demo';
end

% Add source directories to path
addpath(fullfile(fileparts(mfilename('fullpath')), 'src'));
addpath(fullfile(fileparts(mfilename('fullpath')), 'scripts'));

% Load configuration
cfg = config();

% Set random seed
rng(cfg.rng_seed);

switch lower(fig_name)
    case 'demo'
        demo_two_wave(cfg);
    case 'fig6'
        run_script('fig6_7', cfg);
    case 'fig7'
        run_script('fig6_7', cfg);
    case 'fig10'
        run_script('fig10_11', cfg);
    case 'fig11'
        run_script('fig10_11', cfg);
    otherwise
        error('Unknown figure: %s', fig_name);
end
end

function demo_two_wave(cfg)
fprintf('=== SAGE Algorithm Demo: Two-Wave Scenario ===\n');

% Two-wave scenario with different Doppler frequencies
Delta_tau = 200e-9;    % 200 ns delay separation
Delta_phi = 3 * pi/180; % 3 deg azimuth separation
Delta_nu  = 80;        % 80 Hz Doppler separation (nu1=20, nu2=100)

% Resolution limits (Eq. 31)
tau_c = 1 / cfg.B_s;  % ~Tc
phi_c = cfg.lambda / (cfg.M * cfg.d) * 180/pi;  % deg
nu_c  = 1 / (cfg.K * cfg.T);  % Hz

fprintf('Resolution limits:\n');
fprintf('  tau_c = %.2f ns\n', tau_c * 1e9);
fprintf('  phi_c = %.2f deg\n', phi_c);
fprintf('  nu_c  = %.2f Hz\n', nu_c);
fprintf('Wave separation:\n');
fprintf('  Delta_tau = %.2f ns (%.2f * tau_c)\n', Delta_tau * 1e9, Delta_tau / tau_c);
fprintf('  Delta_phi = %.2f deg (%.2f * phi_c)\n', Delta_phi * 180/pi, Delta_phi / (phi_c * pi/180));
fprintf('  Delta_nu  = %.2f Hz (wave1: 20 Hz, wave2: 100 Hz)\n', Delta_nu);

% True parameters
% Wave 1: nu = 20 Hz
alpha_mag = 1;
alpha_phase = 0;
theta1 = [0, 0, 20, alpha_mag * exp(1j * alpha_phase)];

% Wave 2: nu = 100 Hz (alpha_2 = -alpha_1, phase diff pi as in Section V-C)
theta2 = [Delta_tau, Delta_phi, 20 + Delta_nu, alpha_mag * exp(1j * (alpha_phase + pi))];

theta_true = [theta1; theta2];

% Generate received signal
fprintf('\nGenerating received signal...\n');
[y, x_true] = signal_model(cfg, theta_true);

% Compute search grid
fprintf('Computing search grid...\n');
grd = compute_search_grid(cfg);

% Limit grid to reasonable range for demo
grd.tau = 0:cfg.dtau:(2 * Delta_tau + cfg.Tc);
grd.phi = (-10 * pi/180):cfg.dphi:(10 * pi/180);
grd.nu  = (-20):cfg.dnu:(120);

fprintf('Grid sizes: tau=%d, phi=%d, nu=%d\n', ...
    length(grd.tau), length(grd.phi), length(grd.nu));

% Initialize SAGE
fprintf('Initializing via successive cancellation...\n');
theta_init = init_successive_cancellation(y, cfg, grd);

fprintf('Initial estimates:\n');
fprintf('  Wave 1: tau=%.2e, phi=%.2f deg, nu=%.1f, |alpha|=%.3f\n', ...
    theta_init(1,1), theta_init(1,2)*180/pi, theta_init(1,3), abs(theta_init(1,4)));
fprintf('  Wave 2: tau=%.2e, phi=%.2f deg, nu=%.1f, |alpha|=%.3f\n', ...
    theta_init(2,1), theta_init(2,2)*180/pi, theta_init(2,3), abs(theta_init(2,4)));

% Run SAGE
fprintf('\nRunning SAGE algorithm (%d iterations)...\n', cfg.n_iter);
cfg.n_iter = 20;  % Reduce for demo
tic;
[theta_hist, logl_hist] = sage_algorithm(y, cfg, theta_init, grd);
elapsed = toc;
fprintf('Done in %.2f s\n', elapsed);

% Final estimates
theta_final = squeeze(theta_hist(end, :, :));
fprintf('\nFinal estimates:\n');
fprintf('  Wave 1: tau=%.2e, phi=%.2f deg, nu=%.1f, |alpha|=%.3f\n', ...
    theta_final(1,1), theta_final(1,2)*180/pi, theta_final(1,3), abs(theta_final(1,4)));
fprintf('  Wave 2: tau=%.2e, phi=%.2f deg, nu=%.1f, |alpha|=%.3f\n', ...
    theta_final(2,1), theta_final(2,2)*180/pi, theta_final(2,3), abs(theta_final(2,4)));

fprintf('\nTrue parameters:\n');
fprintf('  Wave 1: tau=%.2e, phi=%.2f deg, nu=%.1f, |alpha|=%.3f\n', ...
    theta_true(1,1), theta_true(1,2)*180/pi, theta_true(1,3), abs(theta_true(1,4)));
fprintf('  Wave 2: tau=%.2e, phi=%.2f deg, nu=%.1f, |alpha|=%.3f\n', ...
    theta_true(2,1), theta_true(2,2)*180/pi, theta_true(2,3), abs(theta_true(2,4)));

% Plot log-likelihood convergence
figure;
plot(0:cfg.n_iter, logl_hist, 'b-o', 'LineWidth', 1.5);
xlabel('Iteration cycle');
ylabel('Log-likelihood');
title('SAGE Algorithm Convergence');
grid on;
saveas(gcf, 'results/figures/demo_convergence.png');
savefig(gcf, 'results/figures/demo_convergence.fig');

% Plot parameter convergence
figure;
subplot(2,2,1);
plot(0:cfg.n_iter, abs(squeeze(theta_hist(:,1,4))), 'b-o', ...
     0:cfg.n_iter, abs(squeeze(theta_hist(:,2,4))), 'r-s');
hold on;
plot([0 cfg.n_iter], [abs(theta_true(1,4)) abs(theta_true(1,4))], 'b--');
plot([0 cfg.n_iter], [abs(theta_true(2,4)) abs(theta_true(2,4))], 'r--');
xlabel('Iteration'); ylabel('|alpha|');
legend('Wave 1', 'Wave 2', 'True W1', 'True W2');
grid on; title('Amplitude convergence');

subplot(2,2,2);
plot(0:cfg.n_iter, squeeze(theta_hist(:,1,1))*1e9, 'b-o', ...
     0:cfg.n_iter, squeeze(theta_hist(:,2,1))*1e9, 'r-s');
hold on;
plot([0 cfg.n_iter], [theta_true(1,1) theta_true(1,1)]*1e9, 'b--');
plot([0 cfg.n_iter], [theta_true(2,1) theta_true(2,1)]*1e9, 'r--');
xlabel('Iteration'); ylabel('Delay (ns)');
legend('Wave 1', 'Wave 2', 'True W1', 'True W2');
grid on; title('Delay convergence');

subplot(2,2,3);
plot(0:cfg.n_iter, squeeze(theta_hist(:,1,2))*180/pi, 'b-o', ...
     0:cfg.n_iter, squeeze(theta_hist(:,2,2))*180/pi, 'r-s');
hold on;
plot([0 cfg.n_iter], [theta_true(1,2) theta_true(1,2)]*180/pi, 'b--');
plot([0 cfg.n_iter], [theta_true(2,2) theta_true(2,2)]*180/pi, 'r--');
xlabel('Iteration'); ylabel('Azimuth (deg)');
legend('Wave 1', 'Wave 2', 'True W1', 'True W2');
grid on; title('Azimuth convergence');

subplot(2,2,4);
plot(0:cfg.n_iter, squeeze(theta_hist(:,1,3)), 'b-o', ...
     0:cfg.n_iter, squeeze(theta_hist(:,2,3)), 'r-s');
hold on;
plot([0 cfg.n_iter], [theta_true(1,3) theta_true(1,3)], 'b--');
plot([0 cfg.n_iter], [theta_true(2,3) theta_true(2,3)], 'r--');
xlabel('Iteration'); ylabel('Doppler (Hz)');
legend('Wave 1', 'Wave 2', 'True W1', 'True W2');
grid on; title('Doppler convergence');
saveas(gcf, 'results/figures/demo_params.png');
savefig(gcf, 'results/figures/demo_params.fig');

fprintf('\nDemo complete. Figures saved to results/figures/\n');
end

function run_script(name, cfg)
    run(fullfile('scripts', name));
end
