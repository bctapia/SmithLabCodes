"""
Copyright 2025 Brandon C. Tapia

MIT License
"""

import os


def get_drive(drive='toshiba'):
    """
    Returns the path to the data directory based on the specified drive.
    """

    if drive == 'toshiba':
        data_path = os.path.join(os.sep, 'mnt', 'd')
    elif drive == 'seagate':
        data_path = os.path.join(os.sep, 'mnt', 'e')
    else:
        raise ValueError("Unsupported drive specified.")

    return data_path
