# Structure-Preserving Methods for the Nonlinear Schrodinger Equation

This repository contains a MATLAB project on structure-preserving time
integration for the one-dimensional focusing nonlinear Schrodinger equation
(NLSE)

```math
i\psi_t + \psi_{xx} + \lambda |\psi|^2\psi = 0,
```

with periodic boundary conditions and Fourier spectral spatial discretization.
The code compares standard explicit methods with geometric numerical methods
that preserve important invariants such as mass, energy, or modified energy.

## What Is Implemented

The project implements five NLSE time-stepping methods:

- Explicit Euler, included as an unstable baseline.
- Classical fourth-order Runge-Kutta (RK4).
- Strang split-step Fourier method (SSFM).
- Crank-Nicolson with Delfour-Fortin-Payre-type nonlinear averaging.
- Besse relaxation scheme with modified-energy diagnostics.

It also includes:

- Mass, Hamiltonian energy, momentum, and Besse modified-energy diagnostics.
- Time-convergence experiments.
- Long-time invariant-drift experiments.
- Optional Fourier 2/3-rule de-aliasing.
- A linear Schrodinger harmonic-oscillator demo comparing RK4, Strang splitting,
  Cayley/Crank-Nicolson, and exact spectral propagation.

## Repository Layout

Core setup and diagnostics:

- `nlse_setup.m` builds the spatial grid, wave numbers, initial condition, time
  step count, and optional de-aliasing mask.
- `nlse_diagnostics.m` computes mass, Hamiltonian energy, and momentum.
- `nlse_besse_energy.m` computes the Besse modified energy.
- `run_method.m` runs a stepper and records diagnostic histories.

Time steppers:

- `step_euler.m`
- `step_rk4.m`
- `step_strang.m`
- `step_cn.m`
- `step_besse.m`

Experiment scripts:

- `smoke_test.m` runs a quick invariant-preservation check.
- `compare_methods.m` compares all NLSE methods at one time step.
- `convergence_study.m` estimates empirical time-convergence orders.
- `long_time_drift.m` measures long-time invariant drift.
- `dealiasing_study.m` compares simulations with and without Fourier
  de-aliasing.
- `linear_schroedinger_demo.m` demonstrates structure preservation in a linear
  Schrodinger problem.

Additional project notes and longer writeups are stored in `files/`.


## Main Experiments

Run the method comparison:

```bash
matlab -batch "compare_methods"
```

Run the time-convergence study:

```bash
matlab -batch "convergence_study"
```

Run the long-time drift test:

```bash
matlab -batch "long_time_drift"
```

Run the de-aliasing study:

```bash
matlab -batch "dealiasing_study"
```

Run the linear Schrodinger demo:

```bash
matlab -batch "linear_schroedinger_demo"
```

Most scripts print summary tables and generate diagnostic figures.

## De-Aliasing

The NLSE cubic nonlinearity can transfer energy to unresolved Fourier modes.
This repository includes optional 2/3-rule Fourier filtering.

De-aliasing is disabled by default:

```matlab
P = nlse_setup();
```

To enable it:

```matlab
P = nlse_setup('dealias', true);
```

To change the retained mode fraction:

```matlab
P = nlse_setup('dealias', true, 'dealias_fraction', 2/3);
```

The helper functions `nlse_filter.m` and `nlse_nonlinear.m` provide the common
filtering path used by the time steppers.


## Reproducibility

The scripts are intended to be run directly from the repository root. The setup
function adjusts `tau` slightly so that `Nsteps * tau == T` exactly, which avoids
small final-time mismatches across experiments. The long-time drift test is more expensive because it includes implicit
Crank-Nicolson iterations.
