function [s, t] = sounding_sequence(P, Tc, Ts)
% [s, t] = sounding_sequence(P, Tc, Ts)
% Generate the transmitted burst signal s(t) with rectangular pulses.
%
% Eq. (1) in Fleury et al. (1999)
%
% Inputs:
%   P  - length of PN sequence
%   Tc - chip duration (s)
%   Ts - sampling interval (s)
%
% Outputs:
%   s - P x 1 sampled sounding sequence (normalized to unit power)
%   t - time vector

N_per_chip = round(Tc / Ts);
N_total = P * N_per_chip;

% Generate random binary PN sequence (BPSK)
rng('default');
rng(1);
pn_seq = 2 * randi([0 1], P, 1) - 1;

% Rectangular pulse shaping
s_full = reshape(pn_seq * ones(1, N_per_chip).', N_total, 1);

% Normalize to unit power: ||s||^2 = 1
s = s_full / sqrt(sum(abs(s_full).^2));

t = (0:N_total-1).' * Ts;
end
