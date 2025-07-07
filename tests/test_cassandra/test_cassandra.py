from smithlab import cassandra as css
import os
import numpy as np


def test_get_param():

    prp = os.path.join(os.path.dirname(os.path.abspath(__file__)), "result.out.prp")

    step, prop = css.get_param(prp, "Nmols_2")

    assert step[0] == int(100)
    assert prop[0] == float(0.40000000E+01)
    assert step[-1] == int(2000000)
    assert prop[-1] == float(0.19200000E+03)
    

def test_average_prop():

    prp = os.path.join(os.path.dirname(os.path.abspath(__file__)), "result.out.prp")

    step, prop = css.get_param(prp, "Nmols_2")
    prop_avg, std = css.avg_prop(prop)

    assert np.isclose(prop_avg, 187.00059997)
    assert np.isclose(std, 3.191131143)


def test_loading_lammps_file():

    prp = os.path.join(os.path.dirname(os.path.abspath(__file__)), "result.out.prp")
    lammps_in = os.path.join(os.path.dirname(os.path.abspath(__file__)), "PIM-1_equil_308_1_UA.lmps")
    load, load_std = css.loading(prp, lammps_in)
    
    assert np.isclose(load, 4.06108106869074)
    assert np.isclose(load_std, 0.06930160798746)

def test_loading_given_mass():

    prp = os.path.join(os.path.dirname(os.path.abspath(__file__)), "result.out.prp")
    load, load_std = css.loading(prp, mw=46046.9999999998)

    assert np.isclose(load, 4.06108106869074)
    assert np.isclose(load_std, 0.06930160798746)


def test_concentration():

    prp = os.path.join(os.path.dirname(os.path.abspath(__file__)), "result.out.prp")
    load, load_std = css.loading(prp, mw=46046.9999999998)
    conc, conc_std = css.concentration(load, load_std, 1.0394028)

    assert np.isclose(conc, 94.6117137441344)
    assert np.isclose(conc_std, 1.61453164465678)

# ====================== TESTING WRAPPERS ======================

def test_pressure_calib():
    return


def test_MC_isotherm():
    return


def test_GCMC_isotherm():
    return
