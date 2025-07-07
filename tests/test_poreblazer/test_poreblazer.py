"""
import os
from smithlab import poreblazer as pb


def test_write_xyz():
    lammps_in = os.path.join(os.path.dirname(os.path.abspath(__file__)), "PIM-1_blah.lmps")

    pb.write_xyz(lammps_in)

    return


def test_write_forcefield():
    lammps_in = os.path.join(os.path.dirname(os.path.abspath(__file__)), "PIM-1_blah.lmps")

    pb.write_forcefield(lammps_in)

    return


def test_write_defaults():
    pb.write_defaults()

    return


def test_write_inputs():
    pb.write_input()

    return


def test_write_PB():
    lammps_in = os.path.join(os.path.dirname(os.path.abspath(__file__)), "PIM-1_blah.lmps")

    pb.setup_pb(lammps_in)

    return
"""
