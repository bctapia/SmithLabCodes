import os
from smithlab import poreblazer as pb

"""
def test_write_xyz():
    lammps_in = os.path.join(os.path.dirname(os.path.abspath(__file__)), "PIM-1_blah.lmps")

    assert pb.write_xyz(lammps_in) is None


def test_write_forcefield():
    lammps_in = os.path.join(os.path.dirname(os.path.abspath(__file__)), "PIM-1_blah.lmps")

    assert pb.write_forcefield(lammps_in) is None


def test_write_defaults():

    assert pb.write_defaults() is None


def test_write_inputs():

    assert pb.write_input() is None


def test_setup_pb():
    lammps_in = os.path.join(os.path.dirname(os.path.abspath(__file__)), "PIM-1_blah.lmps")

    assert pb.setup_pb(lammps_in) is None
"""


def test_get_data():
    summary = os.path.join(os.path.dirname(os.path.abspath(__file__)), "summary.dat")
    assert pb.get_data("FV_PO", summary) == 0.29745
    assert pb.get_data("FV_PO", summary, network=True) == 0.23697


def test_get_psd():
    folder = os.path.dirname(os.path.abspath(__file__))
    probe_size_1, psd_1 = pb.get_psd(folder, normalize=False)
    probe_size_2, psd_2 = pb.get_psd(folder, normalize=True)
    probe_size_3, psd_3 = pb.get_psd(folder, cumulative=True, normalize=False)
    assert (probe_size_1[4], psd_1[4]) == (1.125, 2.8245449e-03)
    assert (probe_size_2[4], psd_2[4]) == (1.125, 0.010811137840294721)
    assert (probe_size_3[4], psd_3[4]) == (0.8750000, 0.9998588)


def test_get_probe():
    folder = os.path.dirname(os.path.abspath(__file__))
    assert pb.get_probe(0.01, folder) == 8.375
