# Path Integral Monte Carlo for the 1D Quantum Harmonic Oscillator

This repository contains a Path Integral Monte Carlo (PIMC) simulation in Python for the **1D quantum harmonic oscillator** with several advanced features:

1. **Multiple temperatures** (by varying $\beta$ = $\frac{1}{k_B, T}$).
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
```

## Usage

1. **Clone** or download this repository:
   ```bash
   git clone https://github.com/<your-username>/pimc-quantum-harmonic-oscillator.git
   cd pimc-quantum-harmonic-oscillator

2. **Run** the simulation:
   ```bash
   python pimc_ho.py


This will:
   - Run multiple simulations at different \(\beta\) values in parallel.
   - Display plots for \(\langle x^2 \rangle\) vs. \(\beta\), \(\langle E \rangle\) vs. \(\beta\), and the correlation function at the lowest temperature.
   - Print a summary of results to the terminal.

3. **Customize** parameters:
   - In `pimc_ho.py`, adjust `N`, `n_sweeps`, `step_size`, `thermal_sweeps`, or the list of `betas` as needed (see `worker_run` and the `main()` function).

---

## Interpreting Results

- For the 1D HO at **low temperature** (\(\beta \to \infty\)), we expect:
  - \(\langle x^2 \rangle \approx 0.5\) in the ground state.
  - \(\langle E \rangle \approx 0.5\), the ground-state energy.

- At **higher temperature** (\(\beta\) small), excited states contribute more, increasing \(\langle x^2 \rangle\) and \(\langle E \rangle\).

---

