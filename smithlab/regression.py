import numpy as np

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
