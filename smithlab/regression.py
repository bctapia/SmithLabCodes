import numpy as np


def loss(x, p, mu):
    """ """

    a = x[0]
    b = x[1]

    mu_model = a * np.log(p) + b

    ssr = np.sum((mu - mu_model) ** 2 / mu**2)
    print(ssr)
    return ssr
