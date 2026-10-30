import numpy as np
from scipy.optimize import least_squares
import matplotlib.pyplot as plt

def normalize_decay(x, y, plateau_frac=0.1):
    """
    Normalizes a decay to between [0,1]
    """
    # Estimate plateau as average of last 10% of points
    n_points = len(x)
    n_plateau = max(1, int(n_points * plateau_frac))
    plateau_val = np.median(y[-n_plateau:])

    # Normalize so that y=1 at x=0 and y=0 at plateau
    y_norm = (y - plateau_val) / (y[0] - plateau_val)

    return y_norm


def linearized_kww(t, c2, cut_head=0.01, cut_tail=0.01):
    # Cut off head and tail
    n_points = len(t)
    n_cut_head = int(n_points * cut_head)
    n_cut_tail = int(n_points * cut_tail)
    t_cut = t[n_cut_head:-n_cut_tail]
    c2_cut = c2[n_cut_head:-n_cut_tail]

    return np.log(t_cut), np.log(-np.log(c2_cut))


def loss(x, p, mu):
    """ """

    a = x[0]
    b = x[1]

    mu_model = a * np.log(p) + b

    ssr = np.sum((mu - mu_model) ** 2 / mu**2)
    print(ssr)
    return ssr


def chemical_potential(x, p, mu):
    """Return residuals for mu = a * ln(p / p_0) - b."""
    a, b = x
    p_0 = 1.0  # bar
    mu_model = a * np.log(p / p_0) - b
    return mu_model - mu


def fit_chemical_potential(p, mu, p_uncertainty, initial_guess=(1.0, 0.0), plot=True):
    """
    Fit mu = a * ln(p / p_0) - b using uncertainty in pressure.

    p and p_uncertainty must be in bar.
    p_uncertainty should represent uncertainty in each mean pressure.
    Assumes independent, approximately Gaussian pressure errors and a > 0.
    """
    p = np.asarray(p, dtype=float)
    mu = np.asarray(mu, dtype=float)
    p_uncertainty = np.asarray(p_uncertainty, dtype=float)

    if p.ndim != 1 or p.shape != mu.shape or p.shape != p_uncertainty.shape:
        raise ValueError("Inputs must be one-dimensional arrays of equal length.")

    if len(p) < 2 or np.unique(mu).size < 2:
        raise ValueError("At least two distinct chemical potentials are required.")

    if not all(np.all(np.isfinite(v)) for v in (p, mu, p_uncertainty)):
        raise ValueError("Inputs must contain only finite values.")

    if np.any(p <= 0) or np.any(p_uncertainty <= 0):
        raise ValueError("Pressures and pressure uncertainties must be positive.")

    p_0 = 1.0  # bar

    # Obtain a starting estimate using the original model.
    initial_fit = least_squares(
        lambda x: x[0] * np.log(p / p_0) - x[1] - mu,
        x0=initial_guess,
        bounds=([1e-12, -np.inf], [np.inf, np.inf]),
    )

    if not initial_fit.success:
        raise RuntimeError(f"Initial fit failed: {initial_fit.message}")

    def residuals(x):
        a, b = x
        p_model = p_0 * np.exp((mu + b) / a)
        return (p_model - p) / p_uncertainty

    result = least_squares(
        residuals,
        x0=initial_fit.x,
        bounds=([1e-12, -np.inf], [np.inf, np.inf]),
        x_scale="jac",
    )

    if not result.success:
        raise RuntimeError(f"Fit failed: {result.message}")

    a, b = result.x

    if plot is True:
        p_curve = np.geomspace(p.min(), p.max(), 300)
        mu_curve = a * np.log(p_curve / p_0) - b

        fig, ax = plt.subplots()

        ax.errorbar(
            p,
            mu,
            xerr=p_uncertainty,
            fmt="o",
            capsize=3,
            label="Simulation",
        )

        ax.plot(
            p_curve,
            mu_curve,
            label=f"Fit: a = {a:.4g}, b = {b:.4g}",
        )

        ax.set_xlabel("Pressure (bar)")
        ax.set_ylabel("Chemical potential")
        ax.legend()
        fig.tight_layout()
        plt.show()

    return a, b