# Path Integral Monte Carlo for the 1D Quantum Harmonic Oscillator

This repository contains a Path Integral Monte Carlo (PIMC) simulation in Python for the **1D quantum harmonic oscillator** with several advanced features:

1. **Multiple temperatures** (by varying $\beta$ = $\tfrac{1}{k_{B}\,T}\)$).
2. **Parallelization** with Python's `multiprocessing`.
3. **Measurement of total energy** ($\langle E \rangle$) and ($\langle x^2 \rangle$).
4. **Imaginary-time correlation function** ($C(k) = \langle x_i \, x_{i+k}\rangle$).

## Features

- **Quantum System**: Uses the discretized Euclidean action for the 1D HO with $(\hbar = m = \omega = 1)$.
- **Metropolis Updates**: Local moves on the path, governed by the Euclidean action difference.
- **Scalable**: Easily adjust the number of sweeps, path length (`N`), or range of betas.
- **Parallel**: Each $\beta$ runs in a separate process, speeding up analysis across temperatures.

## Requirements

- **Python 3.7+** (tested up to Python 3.10)
- **numpy**
- **matplotlib**
- *(Optional)* `multiprocessing` is part of standard Python, so no extra install is needed.

Install dependencies (e.g., with `pip`):
```bash
pip install numpy matplotlib
