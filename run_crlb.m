% run_crlb.m - Compute CRLB for two-wave scenario (Figures 6-9)
addpath('src'); addpath('scripts');
cfg = config();

tau_c = 1/cfg.B_s;
phi_c = cfg.lambda/(cfg.M*cfg.d);
nu_c  = 1/(cfg.K*cfg.T);
fprintf('Resolution: tau_c=%.2ens phi_c=%.2fdeg nu_c=%.1fHz\n', ...
    tau_c, phi_c*180/pi, nu_c);

dt_vals = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0] * tau_c;
dp_vals = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0] * phi_c;
dn_vals = [0, 0.5] * nu_c;

alpha_mag = 1;
fprintf('Computing CRLB...\n');
for inu = 1:length(dn_vals)
    dn = dn_vals(inu);
    CRLB_a = zeros(length(dt_vals), length(dp_vals));
    CRLB_p = zeros(length(dt_vals), length(dp_vals));
    for it = 1:length(dt_vals)
        for ip = 1:length(dp_vals)
            dt = dt_vals(it);
            dp = dp_vals(ip);
            theta1 = [0, 0, 0, alpha_mag];
            theta2 = [dt, dp, dn, alpha_mag * exp(1j*pi)];
            crlb = crlb_two_waves(cfg, theta1, theta2);
            CRLB_a(it, ip) = sqrt(crlb.alpha1_mag) / alpha_mag;
            CRLB_p(it, ip) = sqrt(crlb.phi1);
        end
    end
    save(sprintf('results/tables/crlb_fig6_7_nu%d.mat', round(dn/nu_c*10)), ...
        'CRLB_a', 'CRLB_p', 'dt_vals', 'dp_vals');
    fprintf('  nu=%.1f done\n', dn/nu_c);
end

figure;
for inu = 1:length(dn_vals)
    load(sprintf('results/tables/crlb_fig6_7_nu%d.mat', round(dn_vals(inu)/nu_c*10)));
    subplot(2,2,inu*2-1);
    imagesc(dt_vals/tau_c, dp_vals/phi_c, CRLB_a');
    colorbar; set(gca,'YDir','normal');
    xlabel('\Delta\tau/\tau_c'); ylabel('\Delta\phi/\phi_c');
    title(sprintf('CRLB(|a_1|)/|a_1|, \\Delta\\nu=%.1f\\nu_c', dn_vals(inu)/nu_c));

    subplot(2,2,inu*2);
    imagesc(dt_vals/tau_c, dp_vals/phi_c, CRLB_p'*180/pi);
    colorbar; set(gca,'YDir','normal');
    xlabel('\Delta\tau/\tau_c'); ylabel('\Delta\phi/\phi_c');
    title(sprintf('CRLB(\\phi_1) [deg], \\Delta\\nu=%.1f\\nu_c', dn_vals(inu)/nu_c));
end
sgtitle('CRLB for Two-Wave Scenario');
saveas(gcf, 'results/figures/fig6_7_crlb.png');
savefig(gcf, 'results/figures/fig6_7_crlb.fig');
fprintf('Done.\n');
exit
