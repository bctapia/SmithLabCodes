"""smithlab.palsfit"""
import subprocess
from pathlib import Path


class RFCFile:
    """
    Represents an .rfc file with editable settings.

    Parses the 'RESOLUTIONFIT DATA BLOCK 4: RESOLUTION FUNCTION' block:
        header
        components (int)
        res_style (str)
        lifetimes (N floats)
        intensities (N floats)
        res_sigma_style (str)
        sigmas (N floats)
    """

    RESOLUTION_HEADER = "RESOLUTIONFIT DATA BLOCK 4: RESOLUTION FUNCTION"

    def __init__(self, path=None):
        self.path = Path(path) if path else None
        self._raw_lines = []

        self.resolution_components = None
        self.res_lifetime_style = None
        self.res_lifetime = None
        self.res_intensity = None
        self.res_sigma_style = None
        self.res_sigma = None

    @classmethod
    def read(cls, file_in):
        obj = cls(file_in)
        file_in = Path(file_in)

        with open(file_in, "r") as f:
            obj._raw_lines = f.readlines()

        obj.parse_settings()
        return obj

    @staticmethod
    def strip_comment(line):
        # Remove anything after a '#'
        return line.split("#", 1)[0].strip()

    @staticmethod
    def tokenize(lines):
        """
        Turn a list of raw lines into tokens, removing comments and blanks.
        Keeps ordering.
        """
        toks = []
        for raw in lines:
            s = RFCFile.strip_comment(raw)
            if not s:
                continue
            toks.extend(s.split())
        return toks

    def parse_settings(self) -> None:
        """
        Find and parse the resolution block.
        """
    
        # find header line index
        header_idx = None
        for i, raw in enumerate(self._raw_lines):
            if self.RESOLUTION_HEADER in raw:
                header_idx = i
                break

        if header_idx is None:
            raise RuntimeError(f"Could not find header: {self.RESOLUTION_HEADER!r}")

        # Tokenize everything AFTER the header line
        toks = self.tokenize(self._raw_lines[header_idx + 1 :])

        # 1) components
        n = int(float(toks[0]))  # tolerate "2" or "2.0"
        if n <= 0:
            raise RuntimeError(f"Resolution components must be positive; got {n}")

        # 2) style token
        lifetime_style = toks[1]

        pos = 2
        lifetimes = []
        for k in range(n):
            lifetimes.append(float(toks[pos + k]))
        pos += n

        intensities = []
        for k in range(n):
            intensities.append(float(toks[pos + k]))
        pos += n

        sigma_style = toks[pos]
        pos += 1

        sigmas = []
        for k in range(n):
            sigmas.append(float(toks[pos + k]))

        # store
        self.resolution_components = n
        self.res_lifetime_style = lifetime_style
        self.res_lifetime = lifetimes
        self.res_intensity = intensities
        self.res_sigma_style = sigma_style
        self.res_sigma = sigmas


    # ---------------------------
    # Writing
    # ---------------------------
    def write(self, file_out=None):
        """
        Write updated resolution block back to file,
        preserving all other lines exactly.
        """
        required = [
            self.resolution_components,
            self.res_lifetime_style,
            self.res_lifetime,
            self.res_intensity,
            self.res_sigma_style,
            self.res_sigma,
            ]

        # If any required value is None → do nothing
        if any(v is None for v in required):
            print("Resolution block not fully defined — skipping write.")
            return

        file_out = Path(file_out) if file_out else self.path

        header = self.RESOLUTION_HEADER

        # Find block start
        start_idx = None
        for i, line in enumerate(self._raw_lines):
            if header in line:
                start_idx = i
                break

        if start_idx is None:
            raise ValueError("Resolution block header not found.")

        # Build new block text
        n = self.resolution_components

        new_block = []
        new_block.append(f"{header}\n")
        new_block.append(f"{n:10d}\n")
        new_block.append(f"{self.res_lifetime_style}\n")

        # lifetimes
        new_block.append(
            " ".join(f"{x:10.5f}" for x in self.res_lifetime) + "\n"
        )

        # intensities
        new_block.append(
            " ".join(f"{x:10.3f}" for x in self.res_intensity) + "\n"
        )

        new_block.append(f"{self.res_sigma_style}\n")

        # sigmas
        new_block.append(
            " ".join(f"{x:10.5f}" for x in self.res_sigma) + "\n"
        )

        # Replace old block
        # We know structure length = 7 lines total
        end_idx = start_idx + 7

        updated_lines = (
            self._raw_lines[:start_idx]
            + new_block
            + self._raw_lines[end_idx:]
        )

        with open(file_out, "w") as f:
            f.writelines(updated_lines)

    def run(self, out_file=None, exe_path="res19.exe", timeout=None):
        """
        Run res19.exe interactively using this RFC file.

        Parameters
        ----------
        out_file : str or Path, optional
            Output filename. If None, will use same stem with .out.
        exe_path : str
            Path to res19 executable.
        timeout : float, optional
            Timeout in seconds.
        """
        if self.path is None:
            raise ValueError("RFC file has no associated path.")

        rfc_file = Path(self.path)
        if not rfc_file.exists():
            raise FileNotFoundError(f"{rfc_file} not found")

        if out_file is None:
            out_file = rfc_file.with_suffix(".out")

        out_file = Path(out_file)

        process = subprocess.Popen(
            exe_path,
            cwd=rfc_file.parent,   # important for legacy exes
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        # Provide responses to interactive prompts
        input_text = f"{rfc_file.name}\n{out_file.name}\n"

        try:
            stdout, stderr = process.communicate(input_text, timeout=timeout)
        except subprocess.TimeoutExpired:
            process.kill()
            raise RuntimeError("res19 execution timed out.")

        if process.returncode != 0:
            raise RuntimeError(
                f"res19 failed with code {process.returncode}\n"
                f"STDERR:\n{stderr}"
            )

        return out_file, stdout



def pos19():
    """runs pos19.exe
    """
    pass

def res19():
    """runs res19.exe
    """
    pass

# rfc = RFCFile.read("input.rfc")

# print(rfc.temperature)        # Access like attribute
# rfc.temperature = "350"       # Modify
# rfc.settings["pressure"] = 5  # Also works

# rfc.write("modified.rfc")

