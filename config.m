function cfg = config()
% cfg = config()
% Simulation parameters for SAGE algorithm
% Based on Fleury et al. (1999), Section V-D

%% Antenna array
cfg.M          = 10;                % Number of antenna elements (ULA)
cfg.lambda     = 1;                 % Wavelength (normalized)
cfg.d          = cfg.lambda / 2;    % Element spacing
cfg.ula_ref    = 'center';          % Reference point: center of ULA

%% Sounding sequence (PN sequence)
cfg.P          = 255;               % Length of PN sequence
cfg.Tc         = 130e-9;            % Chip duration (130 ns)
cfg.N_c        = 8;                 % Samples per chip (oversampling factor)
cfg.Ts         = cfg.Tc / cfg.N_c;  % Sampling interval (= Tc/8 for oversampling)
cfg.T          = cfg.P * cfg.Tc;    % Burst duration
cfg.pulse      = 'rect';            % Rectangular pulse

%% Observation
cfg.K          = 20;                % Number of snapshots (observation intervals)
cfg.T_obs      = cfg.T;             % Observation duration per snapshot
cfg.T_span     = cfg.K * cfg.T;     % Observation span

%% Quantization precision
cfg.dtau       = cfg.Tc / 8;        % Delay quantization (Tc/8 = 16.25 ns)
cfg.dphi       = 1 * pi/180;        % Azimuth quantization (1 deg)
cfg.dnu        = 5;                 % Doppler quantization (5 Hz)

%% SAGE algorithm parameters
cfg.L          = 2;                 % Number of waves (two-wave scenario)
cfg.n_iter     = 100;               % Number of iteration cycles
cfg.beta       = 1;                 % Data decomposition factor (allocation)

%% Signal-to-noise ratio
cfg.SNR_dB     = 10;                % SNR per wave per antenna branch (dB)
cfg.SNR_lin    = 10^(cfg.SNR_dB/10);

%% Monte Carlo
cfg.Nmc        = 200;               % Number of Monte Carlo runs

%% Random seed
cfg.rng_seed   = 1;

%% Derived parameters
cfg.fs         = 1 / cfg.Ts;        % Sampling frequency
cfg.B_s        = 1 / cfg.Tc;        % Signal bandwidth (approx)
cfg.delay_max  = cfg.P * cfg.Tc;    % Maximum delay search range
cfg.Ns         = round(cfg.P * cfg.Tc / cfg.Ts);  % Samples per burst

%% Generate fixed sounding sequence
rng(42);  % Fixed seed for PN sequence (separate from noise seed)
pn_seq = 2 * randi([0 1], cfg.P, 1) - 1;
N_per_chip = round(cfg.Tc / cfg.Ts);
N_total = cfg.Ns;
s_full = reshape(pn_seq * ones(1, N_per_chip), N_total, 1);
cfg.s = s_full / sqrt(sum(abs(s_full).^2));
cfg.s_power = sum(abs(cfg.s).^2);
