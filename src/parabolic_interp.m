function [peak, peak_val] = parabolic_interp(values, grid)
% [peak, peak_val] = parabolic_interp(values, grid)
% Parabolic interpolation for sub-grid peak refinement (Section IV-A).
%
% Inputs:
%   values - cost function values on uniform grid
%   grid   - grid vector (positions)
%
% Outputs:
%   peak     - refined peak position
%   peak_val - interpolated value at peak

[~, idx] = max(values);

if idx <= 1 || idx >= length(values)
    peak = grid(idx);
    peak_val = values(idx);
    return;
end

y_m1 = values(idx - 1);
y_0  = values(idx);
y_p1 = values(idx + 1);

step = grid(2) - grid(1);

% Parabola: y = a*x^2 + b*x + c  (x in grid units, x=0 at peak)
a = (y_m1 + y_p1 - 2*y_0) / 2;
b = (y_p1 - y_m1) / 2;

if a >= 0
    peak = grid(idx);
    peak_val = y_0;
    return;
end

delta = -b / (2 * a);
delta = max(min(delta, 0.5), -0.5);  % Clamp to neighbor bounds

peak = grid(idx) + delta * step;
peak_val = a * delta^2 + b * delta + y_0;
end
