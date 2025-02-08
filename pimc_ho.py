#!/usr/bin/env python3
"""
pimc_ho.py
----------
Path Integral Monte Carlo (PIMC) for the 1D Quantum Harmonic Oscillator


1) Multiple temperatures (loop over beta = 1/(k_B * T)).
2) Parallelization using Python's multiprocessing (one process per beta).
3) Measurement of total energy (kinetic + potential) as well as <x^2>.
4) Calculation of the imaginary-time correlation function C(k).

Usage:
------
    python pimc_ho.py

This will:
1) Run PIMC simulations for a range of beta (inverse temperature) values.
2) Plot <x^2> vs. beta, plot <E> vs. beta, and plot the correlation function
   for the largest beta.
3) Print a summary of results in the console.
"""

import numpy as np
import matplotlib.pyplot as plt
import multiprocessing


##############################################################################
#                           PIMC CORE FUNCTIONS
##############################################################################

def action_difference(path, i, old_xi, new_xi, dtau):
    """
    Compute the local change in Euclidean action when x_i is replaced by new_xi.

    Parameters
    ----------
    path : np.ndarray
        Current path in imaginary time (shape: (N,))
    i : int
        Index of the slice to be changed
    old_xi, new_xi : float
        Old and proposed new values of x_i
    dtau : float
        Time step (beta / N)

    Returns
    -------
    dS : float
        The change in action (S_new - S_old) associated with flipping x_i.
    """
    N = len(path)
    i_prev = (i - 1) % N
    i_next = (i + 1) % N

    # Old local action
    old_segment = 0.0
    # Kinetic from (i-1 -> i)
    old_segment += 0.5 * (path[i] - path[i_prev])**2 / dtau
    # Kinetic from (i -> i+1)
    old_segment += 0.5 * (path[i_next] - path[i])**2 / dtau
    # Potential at i
    old_segment += 0.5 * dtau * (path[i]**2)

    # New local action
    new_segment = 0.0
    new_segment += 0.5 * (new_xi - path[i_prev])**2 / dtau
    new_segment += 0.5 * (path[i_next] - new_xi)**2 / dtau
    new_segment += 0.5 * dtau * (new_xi**2)

    return new_segment - old_segment


def metropolis_sweep(path, dtau, step_size):
    """
    Perform one Metropolis sweep (N local updates) on the path.

    Parameters
    ----------
    path : np.ndarray
        The current path in imaginary time (shape: (N,))
    dtau : float
        Time step = beta / N
    step_size : float
        Maximum uniform displacement for updates

    Returns
    -------
    path : np.ndarray
        The updated path
    """
    N = len(path)
    for _ in range(N):
        i = np.random.randint(0, N)
        old_xi = path[i]
        new_xi = old_xi + np.random.uniform(-step_size, step_size)

        dS = action_difference(path, i, old_xi, new_xi, dtau)
        if dS < 0.0:
            # Always accept if action decreases
            path[i] = new_xi
        else:
            # Accept with probability exp(-dS)
            if np.random.rand() < np.exp(-dS):
                path[i] = new_xi
    return path


def measure_observables(path, dtau):
    """
    Measure observables from the current path configuration:
    - <x^2>
    - Total energy = Kinetic + Potential
    - Imaginary-time correlation function C(k) = < x_i x_{i+k} >

    Parameters
    ----------
    path : np.ndarray
        The current path (shape: (N,))
    dtau : float
        Time step = beta / N

    Returns
    -------
    dict
        Dictionary containing 'x2', 'energy', and 'correlation'
    """
    N = len(path)

    # 1) <x^2>
    x2 = np.mean(path**2)

    # 2) Estimate the energy using a simple finite-difference approach:
    #    Kinetic ~ 0.5 * <(x_{i+1} - x_i)^2> / dtau^2
    #    Potential ~ 0.5 * <x^2>
    diffs = path - np.roll(path, -1)
    T_est = 0.5 * np.mean(diffs**2) / (dtau**2)
    V_est = 0.5 * x2
    E_est = T_est + V_est

    # 3) Correlation function up to k = N//2
    max_k = N // 2
    corr = np.zeros(max_k)
    for k in range(max_k):
        product = path * np.roll(path, -k)
        corr[k] = np.mean(product)

    return {
        "x2": x2,
        "energy": E_est,
        "correlation": corr
    }


def run_pimc_1d_ho(beta, N=100, n_sweeps=40000, step_size=0.5, thermal_sweeps=8000):
    """
    Run PIMC for the 1D Harmonic Oscillator at a single beta (inverse temperature).

    Parameters
    ----------
    beta : float
        Inverse temperature (1 / (k_B * T)) [k_B = 1 in natural units]
    N : int
        Number of time slices for the path
    n_sweeps : int
        Total number of Metropolis sweeps
    step_size : float
        Proposal step size for uniform spin updates
    thermal_sweeps : int
        Number of sweeps to discard (thermalization period)

    Returns
    -------
    dict
        Contains final averages and histories of:
        - x2_history : list of <x^2> measurements
        - E_history  : list of total energy measurements
        - final_corr : correlation function from the last measurement
        - x2_avg     : final average of <x^2>
        - E_avg      : final average of total energy
    """
    dtau = beta / N
    # Initialize path (random or zeros)
    path = np.random.normal(0.0, 1.0, size=N)

    x2_history = []
    E_history  = []
    final_corr = None

    for sweep in range(n_sweeps):
        path = metropolis_sweep(path, dtau, step_size)

        # Measurements after thermalization
        if sweep >= thermal_sweeps:
            obs = measure_observables(path, dtau)
            x2_history.append(obs["x2"])
            E_history.append(obs["energy"])
            final_corr = obs["correlation"]

    x2_avg = np.mean(x2_history)
    E_avg  = np.mean(E_history)

    return {
        "beta": beta,
        "N": N,
        "x2_history": x2_history,
        "E_history": E_history,
        "final_corr": final_corr,
        "x2_avg": x2_avg,
        "E_avg": E_avg
    }


##############################################################################
#                   PARALLEL DRIVER: MULTIPLE BETAS
##############################################################################

def worker_run(beta):
    """
    Worker function for parallel execution. Adjust parameters as needed.
    """
    # Customize N, n_sweeps, step_size, thermal_sweeps as desired
    N = 100
    n_sweeps = 40000
    step_size = 0.5
    thermal_sweeps = 8000
    result = run_pimc_1d_ho(beta, N, n_sweeps, step_size, thermal_sweeps)
    return result


def main():
    # Example: run at multiple betas => multiple temperatures
    # From beta=1 (T=1) to beta=6 (T~0.1667)
    betas = np.linspace(1.0, 6.0, 6)

    # Run in parallel
    with multiprocessing.Pool(processes=len(betas)) as pool:
        results = pool.map(worker_run, betas)

    # Sort results by beta
    results.sort(key=lambda d: d["beta"])

    # Extract data for plotting
    betas_sorted = [res["beta"] for res in results]
    x2s = [res["x2_avg"] for res in results]
    Es  = [res["E_avg"]  for res in results]

    # Plot <x^2> vs beta
    plt.figure(figsize=(7,5))
    plt.plot(betas_sorted, x2s, 'o--', label=r'$\langle x^2 \rangle$')
    plt.xlabel(r'$\beta$')
    plt.ylabel(r'$\langle x^2 \rangle$')
    plt.title('1D HO: <x^2> vs. Beta')
    plt.grid(True)
    plt.legend()
    plt.show()

    # Plot <E> vs beta
    plt.figure(figsize=(7,5))
    plt.plot(betas_sorted, Es, 'o--', color='red', label=r'$\langle E \rangle$')
    plt.xlabel(r'$\beta$')
    plt.ylabel(r'$\langle E \rangle$')
    plt.title('1D HO: <E> vs. Beta')
    plt.grid(True)
    plt.legend()
    plt.show()

    # Correlation function for the largest beta (lowest T)
    best_result = results[-1]  # largest beta after sorting
    corr = best_result["final_corr"]
    max_k = len(corr)
    plt.figure(figsize=(7,5))
    plt.plot(range(max_k), corr, 'b.-')
    plt.xlabel('Imaginary-time separation k')
    plt.ylabel(r'$C(k) = \langle x_i \, x_{i+k}\rangle$')
    plt.title(f'Correlation Function at beta={best_result["beta"]:.2f}')
    plt.grid(True)
    plt.show()

    # Print summary
    print("=== Summary of Results ===")
    for r in results:
        print(f"Beta = {r['beta']:.2f} | <x^2> = {r['x2_avg']:.4f} | <E> = {r['E_avg']:.4f}")


if __name__ == "__main__":
    main()
