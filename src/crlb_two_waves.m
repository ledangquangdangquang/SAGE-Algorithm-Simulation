function crlb_2 = crlb_two_waves(cfg, theta1, theta2, known_params)
% crlb_2 = crlb_two_waves(cfg, theta1, theta2, known_params)
% Compute CRLB for two-wave scenario by inverting the FIM.
%
% known_params: cell array of parameter names to treat as known
%   e.g., {'nu1','nu2'} means Doppler is known (not estimated)
%
% Uses Eqs. (22)-(25) and Appendix B from Fleury et al. (1999).

if nargin < 4
    known_params = {};
end

M = cfg.M; K = cfg.K; P = cfg.P; Tc = cfg.Tc; Ts = cfg.Ts;
d = cfg.d; lambda = cfg.lambda; Ns = cfg.Ns;
s = cfg.s; s_power = cfg.s_power;

N0 = 2 * abs(theta1(4))^2 * s_power / cfg.SNR_lin;
t = (0:Ns-1) * Ts;
ds_dt = [diff(s); 0] / Ts;

c1 = steering_vector(theta1(2), M, d, lambda);
c2 = steering_vector(theta2(2), M, d, lambda);
m = (0:M-1).' - (M-1)/2;
dc1 = 1j * 2 * pi * d / lambda * cos(theta1(2)) * m .* c1;
dc2 = 1j * 2 * pi * d / lambda * cos(theta2(2)) * m .* c2;

% Parameter labels: tau1, phi1, nu1, Re_a1, Im_a1, tau2, phi2, nu2, Re_a2, Im_a2
param_labels = {'tau1','phi1','nu1','Re_a1','Im_a1','tau2','phi2','nu2','Re_a2','Im_a2'};
n_total = 10;

% Determine which parameters are unknown
unknown = true(1, n_total);
for k = 1:length(known_params)
    idx = find(strcmp(param_labels, known_params{k}));
    if ~isempty(idx)
        unknown(idx) = false;
    end
end
n_unknown = sum(unknown);

% Maps from unknown index to full index
unknown_idx = find(unknown);
full_from_unknown = @(u) unknown_idx(u);

% Derivatives
deriv_cell = cell(n_total, 1);

for k = 1:K
    t_k = t + (k-1)*P*Tc;
    d1 = exp(1j*2*pi*theta1(3)*t_k);
    d2 = exp(1j*2*pi*theta2(3)*t_k);

    s1 = circshift(s, round(theta1(1)/Ts));
    s2 = circshift(s, round(theta2(1)/Ts));
    s1d = circshift(ds_dt, round(theta1(1)/Ts));
    s2d = circshift(ds_dt, round(theta2(1)/Ts));

    x1 = theta1(4) * (c1 * s1.'); x1 = x1 .* repmat(d1, M, 1);
    x2 = theta2(4) * (c2 * s2.'); x2 = x2 .* repmat(d2, M, 1);

    % Wave 1 derivatives
    dtau1 = -theta1(4) * (c1 * s1d.'); dtau1 = dtau1 .* repmat(d1, M, 1);
    dphi1 = theta1(4) * (dc1 * s1.'); dphi1 = dphi1 .* repmat(d1, M, 1);
    dnu1  = theta1(4) * (c1 * s1.'); dnu1 = dnu1 .* repmat(1j*2*pi*t_k.*d1, M, 1);

    % Wave 2 derivatives
    dtau2 = -theta2(4) * (c2 * s2d.'); dtau2 = dtau2 .* repmat(d2, M, 1);
    dphi2 = theta2(4) * (dc2 * s2.'); dphi2 = dphi2 .* repmat(d2, M, 1);
    dnu2  = theta2(4) * (c2 * s2.'); dnu2 = dnu2 .* repmat(1j*2*pi*t_k.*d2, M, 1);

    if k == 1
        for i = 1:n_total
            deriv_cell{i} = zeros(M*Ns*K, 1);
        end
    end

    idx_slice = (k-1)*M*Ns + 1 : k*M*Ns;
    deriv_cell{1}(idx_slice) = dtau1(:);
    deriv_cell{2}(idx_slice) = dphi1(:);
    deriv_cell{3}(idx_slice) = dnu1(:);
    deriv_cell{4}(idx_slice) = x1(:) / theta1(4);
    deriv_cell{5}(idx_slice) = 1j * x1(:) / theta1(4);
    deriv_cell{6}(idx_slice) = dtau2(:);
    deriv_cell{7}(idx_slice) = dphi2(:);
    deriv_cell{8}(idx_slice) = dnu2(:);
    deriv_cell{9}(idx_slice) = x2(:) / theta2(4);
    deriv_cell{10}(idx_slice) = 1j * x2(:) / theta2(4);
end

% Build FIM for unknown parameters
F = zeros(n_unknown, n_unknown);
for ui = 1:n_unknown
    i = full_from_unknown(ui);
    for uj = 1:n_unknown
        j = full_from_unknown(uj);
        F(ui, uj) = (2 / N0) * real(deriv_cell{i}' * deriv_cell{j});
    end
end

% Invert with scaling for numerical stability
scales = [1e-7, 1, 1e4, 1, 1, 1e-7, 1, 1e4, 1, 1];
scales_u = scales(unknown);
S = diag(scales_u);
F_s = S * F * S;
invF_s = inv(F_s + 1e-10 * eye(n_unknown));
invF_u = S * invF_s * S;

% Extract CRLBs for wave 1
crlb_2.tau1 = NaN;
crlb_2.phi1 = NaN;
crlb_2.nu1 = NaN;
crlb_2.alpha1_mag = NaN;

% Map from original parameter indices to unknown indices
orig_to_u = zeros(1, n_total);
for ui = 1:n_unknown
    orig_to_u(unknown_idx(ui)) = ui;
end

% CRLB for tau1
if orig_to_u(1) > 0
    crlb_2.tau1 = invF_u(orig_to_u(1), orig_to_u(1));
end

% CRLB for phi1
if orig_to_u(2) > 0
    crlb_2.phi1 = invF_u(orig_to_u(2), orig_to_u(2));
end

% CRLB for nu1
if orig_to_u(3) > 0
    crlb_2.nu1 = invF_u(orig_to_u(3), orig_to_u(3));
end

% CRLB for |alpha1|
if orig_to_u(4) > 0 && orig_to_u(5) > 0
    ui_re = orig_to_u(4);
    ui_im = orig_to_u(5);
    Re_a = real(theta1(4));
    Im_a = imag(theta1(4));
    alpha_mag = abs(theta1(4));
    crlb_2.alpha1_mag = (Re_a^2 * invF_u(ui_re, ui_re) + ...
        Im_a^2 * invF_u(ui_im, ui_im) + ...
        2*Re_a*Im_a*invF_u(ui_re, ui_im)) / alpha_mag^2;
end

% RMSEE
crlb_2.rmsee_tau1 = sqrt(max(crlb_2.tau1, 0));
crlb_2.rmsee_phi1 = sqrt(max(crlb_2.phi1, 0));
crlb_2.rmsee_nu1 = sqrt(max(crlb_2.nu1, 0));
crlb_2.rmsee_alpha1 = sqrt(max(crlb_2.alpha1_mag, 0));
end
