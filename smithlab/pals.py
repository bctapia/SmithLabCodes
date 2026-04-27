"""
Designed to help analyze data from positron annihilation lifetime spectroscopy (PALS)
"""
import numpy as np


def linearize_dat(dat_in, csv_out):
    """Linearizes the data in a PALSfit3 .dat file and saves it as a .csv file
    """
    counts = []

    with open(dat_in, "r") as f:
        # First line is comment
        comment = f.readline().strip()

        # Read remaining lines
        for line in f:
            if line.strip():  # skip empty lines
                values = line.split()
                counts.extend(int(v) for v in values)

    counts = np.array(counts)

    # Create channel index
    channels = np.arange(len(counts))

    # Save output
    np.savetxt(
        csv_out,
        np.column_stack((channels, counts)),
        delimiter=" ",
        header=comment,
        comments="# ",
        fmt="%d %d"
    )
