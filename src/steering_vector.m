function c_phi = steering_vector(phi, M, d, lambda)
% c_phi = steering_vector(phi, M, d, lambda)
% Steering vector for Uniform Linear Array (ULA).
% Reference point: center of array.
%
% Eq. (1) in Fleury et al. (1999)
%
% Inputs:
%   phi    - incidence azimuth (rad), 0 = broadside
%   M      - number of sensors
%   d      - inter-element spacing (m)
%   lambda - wavelength (m)
%
% Output:
%   c_phi  - M x 1 steering vector

m = (0:M-1).' - (M-1)/2;
c_phi = exp(1j * 2 * pi * d / lambda * m * sin(phi));
end
