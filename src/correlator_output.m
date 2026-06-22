function z = correlator_output(x, tau, phi, nu, cfg)
% z = correlator_output(x, tau, phi, nu, cfg)
% Compute correlator output z(tau, phi, nu; x) from Eq. (12).
%
% z(tau, phi, nu; x) = sum_k sum_m x_m[k] * conj(s(kTs - tau))
%                       * exp(-j*2*pi*nu*kTs) * conj(c_m(phi))
%
% Inputs:
%   x   - M x Nsamp x K array (signal to correlate)
%   tau - test delay (s)
%   phi - test azimuth (rad)
%   nu  - test Doppler (Hz)
%   cfg - configuration struct
%
% Output:
%   z   - scalar complex correlator output

M = cfg.M;
K = cfg.K;
P = cfg.P;
Tc = cfg.Tc;
Ts = cfg.Ts;
d = cfg.d;
lambda = cfg.lambda;
Ns = cfg.Ns;
s = cfg.s;
t = (0:Ns-1).' * Ts;

% Steering vector for test azimuth
c_phi = steering_vector(phi, M, d, lambda);

% Delayed sounding sequence
tau_samp = round(tau / Ts);
s_delayed = circshift(s, tau_samp);

% Compute z
z = 0;
for k = 1:K
    % Time base
    t_center = (k - 1) * P * Tc;
    t_k = t + t_center;

    % Doppler compensation
    doppler = exp(-1j * 2 * pi * nu * t_k);

    % Compute: (c_phi^H * x[:,:,k]) · conj(s_delayed) · exp(-j*2*pi*nu*t_k)
    beamformed = c_phi' * x(:, :, k);  % 1 x Ns
    template = conj(s_delayed.') .* doppler.';
    z_k = sum(beamformed .* template);
    z = z + z_k;
end
end
