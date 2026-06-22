# SAGE Algorithm Simulation

MATLAB implementation of the **SAGE (Space-Alternating Generalized Expectation-Maximization)** algorithm from Fleury et al. (1999) — joint estimation of delay, azimuth, Doppler frequency, and complex amplitude of multipath waves in mobile radio channels.

## Structure

```
SAGE/
├── main.m          # Entry point: demo + figure dispatch
├── config.m        # Simulation parameters
├── src/            # Core functions
│   ├── signal_model.m                # Received signal generation (Eq. 2)
│   ├── steering_vector.m             # ULA steering vector (Eq. 1)
│   ├── correlator_output.m           # Correlation z(τ,φ,ν; x) (Eq. 12)
│   ├── mle_single_wave.m             # M-step: coordinate-wise MLE + parabolic interp
│   ├── expectation_step.m            # E-step (Eq. 13)
│   ├── sage_algorithm.m              # Main SAGE loop
│   ├── init_successive_cancellation.m # Initialization (Section IV-B)
│   ├── parabolic_interp.m            # Parabolic sub-grid refinement
│   ├── compute_search_grid.m         # τ, φ, ν search grids
│   ├── log_likelihood.m              # Log-likelihood (Eq. 4)
│   ├── crlb_single_wave.m            # Single-wave CRLB (Eqs. 27-30)
│   └── crlb_two_waves.m              # Two-wave CRLB
├── scripts/        # Figure reproduction
│   ├── fig6_7.m
│   └── fig10_11.m
├── results/        # Output
│   ├── figures/
│   └── tables/
└── docs/           # Documentation (AI-generated)
```

## Usage

```matlab
main           % Two-wave demo
main('demo')   % Same as above
main('fig6')   % CRLB Figures 6-7
main('fig10')  % Convergence Figures 10-11
```

## Key Parameters (config.m)

| Parameter | Value | Description |
|-----------|-------|-------------|
| M | 10 | ULA elements |
| P | 255 | PN sequence length |
| Tc | 130 ns | Chip duration |
| N_c | 8 | Samples per chip (oversampling) |
| Ts | 16.25 ns | Sampling interval |
| K | 20 | Snapshots |
| SNR | 10 dB | SNR per wave per antenna |
| n_iter | 100 | SAGE iterations |

## Demo results

With Δτ = 200ns (>τ_c), Δφ = 3° (<φ_c/2):

```
Wave 1: τ=0.0ns, φ=0°,   ν=16.3Hz, |α|=1.042  (true: 0, 0°, 20Hz, 1.0)
Wave 2: τ=195ns, φ=2.9°, ν=116.9Hz,|α|=0.981  (true: 200ns, 3°, 100Hz, 1.0)
```

SAGE resolves both waves even though azimuth separation (3°) is below the classical resolution φ_c/2 = 5.73°. Waves are separable through delay (Δτ > τ_c) per Eq. (32).

## Reference

Fleury, B. H., Tschudin, M., Heddergott, R., Dahlhaus, D., & Pedersen, K. I. (1999). Channel parameter estimation in mobile radio environments using the SAGE algorithm. *IEEE Journal on Selected Areas in Communications*, 17(3), 434–450.
