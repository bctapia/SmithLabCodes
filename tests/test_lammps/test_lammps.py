from smithlab import lammps
import os
import numpy as np


def test_reformat():

    lammps_in = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "PIM-1_equil_unformat_308_1_UA.lmps"
    )
    lammps_out = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "PIM-1_equil_308_1_UA.lmps"
    )
    lammps_ref = os.path.join(os.path.dirname(os.path.abspath(__file__)), "PIM-1_MON.lmps")

    assert lammps.reformat(lammps_in, lammps_out, lammps_ref) is None


def test_reorder_sections():

    lammps_in = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "PIM-1_equil_unformat_308_1_UA.lmps"
    )
    lammps_out = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "PIM-1_equil_308_1_UA.lmps"
    )

    assert lammps.reorder_sections(lammps_in, lammps_out) is None


def test_setup_lammps():

    lammps_in = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "PIM-1_equil_unformat_308_1_UA.lmps"
    )
    lammps_out = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "PIM-1_equil_308_1_UA.lmps"
    )
    lammps_ref = os.path.join(os.path.dirname(os.path.abspath(__file__)), "PIM-1_MON.lmps")

    assert lammps.setup_lammps(lammps_in, lammps_out, lammps_ref) is None


def test_polym_stats():

    lammps_in = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "PIM-1_equil_308_1_UA.lmps"
    )
    molecule, count = lammps.polym_stats(lammps_in)

    assert np.array_equal(molecule, np.array([1, 2, 3, 4]))
    assert np.array_equal(count, np.array([1575, 70, 945, 910]))


def test_get_mw():

    lammps_in = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "PIM-1_equil_308_1_UA.lmps"
    )

    assert lammps.get_mw(lammps_in) == 46047.0
