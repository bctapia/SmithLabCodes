"""
Copyright 2025 Brandon C. Tapia

MIT License
"""

from pathlib import Path
import platform


def get_drive(drive="d"):
    """
    Returns a Path object to the specified drive.

    Returns:
    - pathlib.Path object
    """

    if platform.system() == "Windows":
        path = Path(f"{drive.upper()}:\\")
    elif platform.system() == "darwin":  # macOS
        path = Path(f"/Volumes/{drive}") 
    else:  # assume Unix-like system (including WSL)
        path = Path("/mnt") / drive

    return path


def get_wsl_home(distro="Ubuntu", user="btapia"):
    return Path(f"\\\\wsl$\\{distro}\\home\\{user}")
