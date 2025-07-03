from smithlab import cassandra as css
import os


def test_get_param():

    prp = os.path.join(os.path.dirname(os.path.abspath(__file__)), "result.out.prp")

    step, prop = css.get_param(prp, "Nmols_1")


def test_average_prop():

    prp = os.path.join(os.path.dirname(os.path.abspath(__file__)), "result.out.prp")

    step, prop = css.get_param(prp, "Nmols_1")
    prop_avg, std = css.avg_prop(step)


def test_pressure_calib():
    return


def test_MC_isotherm():
    return


def test_GCMC_isotherm():
    return
