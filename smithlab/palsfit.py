"""smithlab.palsfit"""
import subprocess
import re
from pathlib import Path


class RFCFile:
    """Represents a .rfc file with editable settings.
    """

    RESOLUTION_HEADER_1 = "RESOLUTIONFIT DATA BLOCK 1: OUTPUT OPTIONS"
    RESOLUTION_HEADER_2 = "RESOLUTIONFIT DATA BLOCK 2: SPECTRUM"
    RESOLUTION_HEADER_3 = "RESOLUTIONFIT DATA BLOCK 3: CHANNEL RANGES. TIME SCALE. TIME-ZERO."
    RESOLUTION_HEADER_4 = "RESOLUTIONFIT DATA BLOCK 4: RESOLUTION FUNCTION"
    RESOLUTION_HEADER_5 = "RESOLUTIONFIT DATA BLOCK 5: LIFETIMES AND INTENSITY CONSTRAINTS"
    RESOLUTION_HEADER_6 = "RESOLUTIONFIT DATA BLOCK 6: BACKGROUND CONSTRAINTS"

    def __init__(self, path=None):
        self.path = Path(path) if path else None
        self._raw_lines = []

        # Block 1 data (Output options)
        # Represented as an block of 0 (False) and 1 (True) (e.g., 0101)
        # Default: 0000
        self.echo = None
        self.iteration = None
        self.resid_plot = None
        self.corr_matrix = None

        # There is also a hidden log-normal fineness that can be applied if desired
        # Default: 32 (e.g., 0000 32 is the same as 0000 which is the same as 0000 0)
        self.log_normal_fineness = None
        
        # Block 2 data (Spectrum)
        self.num_channels = None  # Number of channels in the spectrum
        self.format = None  # Formatting of spectrum expressed in FORMAT style of FORTRAN
        self.spectrum_path = None  # File path of the spectrum
        self.spectrum_label = None  # Spectrum label (header row in the spectrum)
        self.inspec = None  # Whether the spectrum is an intrinsic part of the spectrum (1) or not (0)

        # Block 3 data (CHANNEL RANGES. TIME SCALE. TIME-ZERO)
        self.area = None  # Integer array corresponding to the min and max area channels [area_min, area_max]
        self.fit = None  # Integer array corresponding to the min and max fit channels [fit_min, fit_max]
        self.timescale = None  # ns/ch
        self.timezero = None  # in channel number

        # Block 4 data (RESOLUTION FUNCTION)
        self.res_components = None  # Number of lifetime components in the resolution func.
        self.res_lt_constraint = None  # Whether components are guessed or fixed (e.g., "GFF")
        self.res_fwhm = None  # Array of FWHM (ns)
        self.res_intensity = None  # Intensity of lifetimes (%)
        self.res_shift_constraint = None  # Whether shifts are guessed or fixed (e.g., "GFF")
        self.res_shift = None  # Array of shifts (peak displacements) (ns)

        # Block 5 data (LIFETIMES AND INTENSITY CONSTRAINTS)
        self.lt_components = None  # Number of lifetime components in the material
        self.lt_constraint = None  # Whether components are guessed or fixed (e.g., "GFF")
        self.lt = None  # Lifetimes (ns)
        self.int_constraint = None # If m=0, relative intensities fixed, if m>0, TODO
        self.int_constraint_info = None  # TODO

        # Block 6 data (BACKGROUND CONSTRAINTS)
        self.bg_constraint = None  # 0=No constraint, 1=fit area between channels in self.bg_channels, 2=fit area specified in self.bg_fixed
        self.bg_channels = None
        self.bg_fixed = None


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
        """Turn a list of raw lines into tokens, removing comments and blanks.
        Keeps ordering.
        """
        toks = []
        for raw in lines:
            s = RFCFile.strip_comment(raw)
            if not s:
                continue
            toks.extend(s.split())
        return toks

    def parse_settings(self):
        """Find and parse the resolution block.
        """

        # find header line index
        header_idx_1 = None
        header_idx_2 = None
        header_idx_3 = None
        header_idx_4 = None
        header_idx_5 = None
        header_idx_6 = None

        for i, raw in enumerate(self._raw_lines):
            if self.RESOLUTION_HEADER_1 in raw:
                header_idx_1 = i
            elif self.RESOLUTION_HEADER_2 in raw:
                header_idx_2 = i
            elif self.RESOLUTION_HEADER_3 in raw:
                header_idx_3 = i
            elif self.RESOLUTION_HEADER_4 in raw:
                header_idx_4 = i
            elif self.RESOLUTION_HEADER_5 in raw:
                header_idx_5 = i
            elif self.RESOLUTION_HEADER_6 in raw:
                header_idx_6 = i

        #==========================BLOCK 1============================
        toks = self.tokenize(self._raw_lines[header_idx_1 + 1:header_idx_2])

        self.echo = int(toks[0][0])
        self.iteration = int(toks[0][1])
        self.resid_plot = int(toks[0][2])
        self.corr_matrix = int(toks[0][3])
        self.log_normal_fineness = int(toks[0][4:].strip())

        #==========================BLOCK 2============================
        toks = self.tokenize(self._raw_lines[header_idx_3 + 1:header_idx_4])

        self.num_channels = toks[0]
        self.format = toks[1]
        self.spectrum_path = toks[2] 
        self.spectrum_label = toks[3]
        self.inspec = toks[4]

        #==========================BLOCK 3============================
        toks = self.tokenize(self._raw_lines[header_idx_3 + 1:header_idx_4])
        
        self.area = [toks[0], toks[1]]
        self.fit = [toks[2], toks[3]]
        self.timescale = toks[4]
        self.timezero = toks[5]

        #==========================BLOCK 4============================
        toks = self.tokenize(self._raw_lines[header_idx_4 + 1:header_idx_5])

        self.res_components = int(float(toks[0]))
        self.res_lt_constraint = toks[1]
        
        self.res_fwhm = []
        pos = 2
        for k in range(self.res_components):
            self.res_fwhm.append(float(toks[pos + k]))
        pos += self.res_components

        self.res_intensity = []
        for k in range(self.res_intensity):
            self.res_intensity.append(float(toks[pos + k]))
        pos += self.res_components

        self.res_shift_constraint = toks[pos]
        pos += 1

        self.res_shift = []
        for k in range(self.res_components):
            self.res_shift.append(float(toks[pos + k]))

        #==========================BLOCK 5============================
        # Tokenize everything AFTER the header line
        toks = self.tokenize(self._raw_lines[header_idx_5 + 1:header_idx_6])

        # 1) components
        n = int(float(toks[0]))  # tolerate "2" or "2.0"
        if n <= 0:
            raise RuntimeError(f"Lifetime components must be positive; got {n}")
        
        # 2) style token
        lifetime_style = toks[1]
        
        pos = 2
        lifetimes = []
        for k in range(n):
            lifetimes.append(float(toks[pos + k]))

        # store
        self.lifetime_components = n
        self.lifetime_style = lifetime_style.strip()
        self.lifetime = lifetimes

        #==========================BLOCK 6============================

    def write(self, file_out=None):
        """
        Write updated resolution block back to file, preserving all other lines exactly.
        """

        file_out = Path(file_out) if file_out else self.path
        self.path = Path(file_out)

        header_1 = self.RESOLUTION_HEADER_1
        header_2 = self.RESOLUTION_HEADER_2
        header_3 = self.RESOLUTION_HEADER_3
        header_4 = self.RESOLUTION_HEADER_4
        header_5 = self.RESOLUTION_HEADER_5
        header_6 = self.RESOLUTION_HEADER_6

        # Find block start
        start_idx_1 = None
        start_idx_2 = None
        start_idx_3 = None
        start_idx_4 = None
        start_idx_5 = None
        start_idx_6 = None
        for i, line in enumerate(self._raw_lines):
            if header_1 in line:
                start_idx_1 = i
            elif header_2 in line:
                start_idx_2 = i
            elif header_3 in line:
                start_idx_3 = i
            elif header_4 in line:
                start_idx_4 = i
            elif header_5 in line:
                start_idx_5 = i
            elif header_6 in line:
                start_idx_6 = i

        #==========================BLOCK 1============================
        new_block = []
        new_block.append(f"{header_1}\n")
        new_block.append(f"{"".join[self.echo, self.iteration, self.resid_plot, self.corr_matrix]}")
        if self.log_normal_fineness:
            new_block.append(f" {self.log_normal_fineness}")
        new_block.append("\n")

        #==========================BLOCK 2============================
        new_block = []
        new_block.append(f"{header_2}\n")
        new_block.append(f"{self.num_channels:>10}\n")
        new_block.append(f"{self.format}\n")
        new_block.append(f"{self.spectrum_path}\n")
        new_block.append(f"{self.spectrum_label}\n")
        new_block.append(f"{self.inspec:>2}\n")
        if int(self.inspec) == 1:
            new_block.append(f"{self.spectrum_label}\n")

        #==========================BLOCK 3============================
        new_block = []
        new_block.append(f"{header_3}\n")
        new_block.append(f"{self.area[0]:>10d}\n")
        new_block.append(f"{self.area[1]:>10d}\n")
        new_block.append(f"{self.fit[0]:>10d}\n")
        new_block.append(f"{self.fit[1]:>10d}\n")
        new_block.append(f"{self.timescale:>10f}\n")
        new_block.append(f"{self.timezero:>10.3f}\n")


        #==========================BLOCK 4============================
        new_block = []
        new_block.append(f"{header_4}\n")
        new_block.append(f"{self.res_components:>10}\n")
        new_block.append(f"{self.res_lt_constraint}")
        new_block.append(" ".join(f"{x:>10.5f}" for x in self.res_fwhm) + "\n")
        new_block.append(" ".join(f"{x:>10.3f}" for x in self.res_intensity) + "\n")
        new_block.append(f"{self.res_shift_constraint}\n")
        new_block.append(" ".join(f"{x:>10.5f}" for x in self.res_shift) + "\n")

        # Replace old block
        # We know structure length = 7 lines total
        #end_idx = start_idx_4 + 7

        #updated_lines = (self._raw_lines[:start_idx_4] + new_block + self._raw_lines[end_idx:])

        #==========================BLOCK 5============================
        new_block = []
        new_block.append(f"{header_5}\n")
        new_block.append(f"{self.lt_components:>10d}\n" )
        new_block.append(f"{self.lt_constraint}\n")
        new_block.append(" ".join(f"{x:>10.5f}" for x in self.lt) + "\n")
        new_block.append(f"{self.int_constraint:>10d}\n")
        if self.int_constraint > 0:
            print("TODO")
        if self.int_constraint < 0:
            print("TODO")
        

        # Replace old block
        # We know structure length = 4 lines total
        #end_idx = start_idx_5 + 4
        #updated_lines = (updated_lines[:start_idx_5] + new_block + updated_lines[end_idx:])

        #with open(file_out, "w") as f:
        #    f.writelines(updated_lines)

        #==========================BLOCK 6============================
        new_block = []
        new_block.append(f"{header_6}\n")
        new_block.append(f"{self.bg_constraint:>10d}\n")
        if int(self.bg_constraint) == 1:
            print("TODO")
        elif int(self.bg_constraint) == 2:
            print("TODO")


    def run(self, exe_path, out_file=None, timeout=None):

        if self.path is None:
            raise ValueError("RFCFile.path is None.")

        exe_path = Path(exe_path)
        rfc_path = Path(self.path)

        if not exe_path.exists():
            raise FileNotFoundError(f"Executable not found: {exe_path}")

        if not rfc_path.exists():
            raise FileNotFoundError(f"RFC not found: {rfc_path}")

        # Default output name in SAME directory as RFC
        if out_file is None:
            out_file = rfc_path.with_suffix(".out")

        out_file = Path(out_file)

        proc = subprocess.Popen([str(exe_path)], cwd=str(rfc_path.parent), stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,)

        # Only send filenames (not full paths)
        stdin_text = f"{rfc_path.name}\n{out_file.name}\n"

        stdout, stderr = proc.communicate(stdin_text, timeout=timeout)

        if proc.returncode != 0:
            raise RuntimeError(f"res19 failed (code {proc.returncode})\nSTDOUT:\n{stdout}\nSTDERR:\n{stderr}")

        produced = rfc_path.parent / out_file.name

        if not produced.exists():
            raise RuntimeError(f"Output file not created: {produced}")

        return produced


class PFCFile:
    """Represents a .pfc file with editable settings
    """

    POSITRON_HEADER_1 = "POSITRONFIT DATA BLOCK 1: OUTPUT OPTIONS"
    POSITRON_HEADER_2 = "POSITRONFIT DATA BLOCK 2: SPECTRUM"
    POSITRON_HEADER_3 = "POSITRONFIT DATA BLOCK 3: CHANNEL RANGES. TIME SCALE. TIME-ZERO."
    POSITRON_HEADER_4 = "POSITRONFIT DATA BLOCK 4: RESOLUTION FUNCTION"
    POSITRON_HEADER_5 = "POSITRONFIT DATA BLOCK 5: LIFETIMES AND INTENSITY CONSTRAINTS"
    POSITRON_HEADER_6 = "POSITRONFIT DATA BLOCK 6: BACKGROUND CONSTRAINTS"
    POSITRON_HEADER_8 = "POSITRONFIT DATA BLOCK 8: SOURCE CORRECTION"

    def __init__(self, path=None):
        self.path = Path(path) if path else None
        self._raw_lines = []

        # Block 1 data (Output options)
        # Represented as an block of 0 (False) and 1 (True) (e.g., 0101)
        # Default: 0000
        self.echo = None
        self.iteration = None
        self.resid_plot = None
        self.corr_matrix = None

        # There is also a hidden log-normal fineness that can be applied if desired
        # Default: 32 (e.g., 0000 32 is the same as 0000 which is the same as 0000 0)
        self.log_normal_fineness = None

        # Block 2 data (Spectrum)
        self.num_channels = None  # Number of channels in the spectrum
        self.format = None  # Formatting of spectrum expressed in FORMAT style of FORTRAN
        self.spectrum_path = None  # File path of the spectrum
        self.spectrum_label = None  # Spectrum label (header row in the spectrum)
        self.inspec = None  # Whether the spectrum is an intrinsic part of the spectrum (1) or not (0)

        # Block 3 data (CHANNEL RANGES. TIME SCALE. TIME-ZERO)
        self.area = None  # Integer array corresponding to the min and max area channels [area_min, area_max]
        self.fit = None  # Integer array corresponding to the min and max fit channels [fit_min, fit_max]
        self.timescale = None  # ns/ch
        self.timezero_constraint = None
        self.timezero = None  # in channel number

        # Block 4 data (RESOLUTION FUNCTION)
        self.res_components = None # Number of lifetime components in the resolution func.
        self.res_fwhm = None # Array of FWHM (ns)
        self.res_intensity = None # Intensity of lifetimes (%)
        self.res_shift = None # Array of shifts (peak displacements) (ns)

        # Block 5 data (LIFETIMES AND INTENSITY CONSTRAINTS)
        self.lt_components = None  # Number of lifetime components in the material
        self.lt_constraint = None  # Whether components are guessed or fixed (e.g., "GFF")
        self.lt = None  # Lifetimes (ns)
        self.ln_constraint = None  # Whether log-normal components are guessed or fixed (e.g., "GFF")
        self.ln_broadening = None  # (ns)
        self.int_constraint = None # If m=0, relative intensities fixed, if m>0, TODO
        self.int_constraint_info = None  # TODO
        self.lt_components_2 = None  # Number of lifetime components in the material
        self.lt_constraint_2 = None  # Whether components are guessed or fixed (e.g., "GFF")
        self.lt_2 = None  # Lifetimes (ns)
        self.ln_constraint_2 = None  # Whether log-normal components are guessed or fixed (e.g., "GFF")
        self.ln_broadening_2 = None  # (ns)
        self.int_constraint_2 = None # If m=0, relative intensities fixed, if m>0, TODO
        self.int_constraint_info_2 = None  # TODO

        # Block 6 data (BACKGROUND CONSTRAINTS)
        self.bg_constraint = None  # 0=No constraint, 1=fit area between channels in self.bg_channels, 2=fit area specified in self.bg_fixed
        self.bg_channels = None
        self.bg_fixed = None

        # Block 7 data (AREA CONSTRAINTS)
        self.area_constraint = None # 0=No constraint, 1=fit area between channels in self.area_channels, 2=fit area specified in self.area_fixed
        self.area_channels = None
        self.area_fixed = None

        # Block 8 data (SOURCE CORRECTION)
        self.source_components = None  # Number of source correction components
        self.source_lt = None  # Lifetimes (ns)
        self.source_broadening = None  # (ns)
        self.source_int = None  # Relative intensities of the source correction components (%)
        self.source_total = None  # Percentage of positrons that annihilate in the source (%)
        self.new_cycle_input = None  # Whether to start source correction convergence from values obtained from no source correction convergence (0) or use new specified params. (1) and whether to also change t0 (2)
        self.source_components_2 = None  # Number of source correction components
        self.source_lt_2 = None  # Lifetimes (ns)
        self.source_broadening_2 = None  # (ns)
        self.source_int_2 = None  # Relative intensities of the source correction components (%)
        self.source_total_2 = None  # Percentage of positrons that annihilate in the source (%)
        self.timezero_constraint_2 = None
        self.timezero_2 = None

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
        """Turn a list of raw lines into tokens, removing comments and blanks.
        Keeps ordering.
        """
        toks = []
        for raw in lines:
            s = PFCFile.strip_comment(raw)
            if not s:
                continue
            toks.extend(s.split())
        return toks

    def parse_settings(self):
        """Find and parse the resolution block.
        """

        header_idx_1 = None
        header_idx_2 = None
        header_idx_3 = None
        header_idx_4 = None
        header_idx_5 = None
        header_idx_6 = None
        header_idx_7 = None
        header_idx_8 = None
        for i, raw in enumerate(self._raw_lines):
            if self.POSITRON_HEADER_1 in raw:
                header_idx_1 = i
            elif self.POSITRON_HEADER_2 in raw:
                header_idx_2 = i
            elif self.POSITRON_HEADER_3 in raw:
                header_idx_3 = i
            elif self.POSITRON_HEADER_4 in raw:
                header_idx_4 = i
            elif self.POSITRON_HEADER_5 in raw:
                header_idx_5 = i
            elif self.POSITRON_HEADER_6 in raw:
                header_idx_6 = i
            elif self.POSITRON_HEADER_7 in raw:
                header_idx_7 = i
            elif self.POSITRON_HEADER_8 in raw:
                header_idx_8 = i

        #==========================BLOCK 1============================
        toks = self.tokenize(self._raw_lines[header_idx_1 + 1:header_idx_2])

        self.echo = int(toks[0][0])
        self.iteration = int(toks[0][1])
        self.resid_plot = int(toks[0][2])
        self.corr_matrix = int(toks[0][3])
        self.log_normal_fineness = int(toks[0][4:].strip())

        #==========================BLOCK 2============================
        toks = self.tokenize(self._raw_lines[header_idx_3 + 1:header_idx_4])

        self.num_channels = toks[0]
        self.format = toks[1]
        self.spectrum_path = toks[2] 
        self.spectrum_label = toks[3]
        self.inspec = toks[4]

        #==========================BLOCK 3============================
        toks = self.tokenize(self._raw_lines[header_idx_3 + 1:header_idx_4])
        
        self.area = [toks[0], toks[1]]
        self.fit = [toks[2], toks[3]]
        self.timescale = toks[4]
        self.timezero_constraint = toks[5]
        self.timezero = toks[6]

        #==========================BLOCK 4============================
        # Tokenize everything AFTER the header line
        toks = self.tokenize(self._raw_lines[header_idx_4 + 1 :header_idx_5])

        self.res_components = int(float(toks[0]))

        pos = 1
        self.res_fwhm = []
        for k in range(self.res_components):
            self.res_fwhm.append(float(toks[pos + k]))
        pos += self.res_components

        self.res_intensity = []
        for k in range(self.res_components):
            self.res_intensity.append(float(toks[pos + k]))
        pos += self.res_components

        self.res_shift = []
        for k in range(self.res_components):
            self.res_shift.append(float(toks[pos + k]))

        #==========================BLOCK 5============================
        toks = self.tokenize(self._raw_lines[header_idx_5 + 1 :])

        self.lt_components = int(float(toks[0]))
        self.lt_constraint = toks[1]
        
        pos = 2
        self.lt = []
        for k in range(self.lt_components):
            self.lt.append(float(toks[pos + k]))
        pos += self.lt_components

        ln_constraint = toks[pos]

        pos += 1
        self.ln_broadening = []
        for k in range(self.lt_components):
            self.ln_broadening.append(float(toks[pos + k]))

        # TODO: figure out potential second iter needed

        #==========================BLOCK 6============================

        #==========================BLOCK 7============================

        #==========================BLOCK 8============================
        # Tokenize everything AFTER the header line
        toks = self.tokenize(self._raw_lines[header_idx_8 + 1 :])

        # 1) components
        n = int(float(toks[0])) # tolerate "2" or "2.0"
        self.source_components = n

        if n > 0:
            pos = 1
            lifetimes = []
            for k in range(n):
                lifetimes.append(float(toks[pos + k]))
            pos += n

            sigmas = []
            for k in range(n):
                sigmas.append(float(toks[pos + k]))
            pos += n

            intensities = []
            for k in range(n):
                intensities.append(float(toks[pos + k]))
            pos += n

            total = float(toks[pos])

            # store
            self.source_lifetime = lifetimes
            self.source_sigma = sigmas
            self.source_intensity = intensities
            self.source_total = total

    # ---------------------------
    # Writing
    # ---------------------------
    def write(self, file_out=None):
        """Write updated resolution block back to file, preserving all other lines exactly.
        """

        file_out = Path(file_out) if file_out else self.path
        self.path = Path(file_out)
        header_4 = self.POSITRON_HEADER_4
        header_5 = self.POSITRON_HEADER_5
        header_6 = self.POSITRON_HEADER_6
        header_8 = self.POSITRON_HEADER_8

        # FInd block start
        start_idx_4 = None
        start_idx_5 = None
        start_idx_6 = None
        start_idx_8 = None
        for i, line in enumerate(self._raw_lines):
            if header_4 in line:
                start_idx_4 = i
            elif header_5 in line:
                start_idx_5 = i
            elif header_6 in line:
                start_idx_6 = i
            elif header_8 in line:
                start_idx_8 = i
        
        # ===================BLOCK 4======================
        new_block = []
        new_block.append(f"{header_4}\n")
        new_block.append(f"{self.resolution_components:10d}\n")

        # lifetimes
        new_block.append(" ".join(f"{x:10.5f}" for x in self.res_lifetime) + "\n")

        # intensities
        new_block.append(" ".join(f"{x:10.3f}" for x in self.res_intensity) + "\n")

        # sigmas
        new_block.append(" ".join(f"{x:10.5f}" for x in self.res_sigma) + "\n")

        # Replace old block
        # We know structure length = 5 lines total
        end_idx = start_idx_4 + 5

        updated_lines = (self._raw_lines[:start_idx_4] + new_block + self._raw_lines[end_idx:])


        # ===================BLOCK 5======================
        new_block = []
        new_block.append(f"{header_5}\n")
        new_block.append(f"{self.lt_components:10d}\n")
        new_block.append(f"{self.lt_style}\n")
        new_block.append(" ".join(f"{x:10.5f}" for x in self.lt_lifetime) + "\n")
        new_block.append(f"{self.ln_style}\n")
        new_block.append(" ".join(f"{x:10.5f}" for x in self.ln_sigma) + "\n")
        
        end_idx = start_idx_6 - 1  # TODO UPDATE THIS ONCE WE ADD OTHER CONSTRAINTS

        updated_lines = (updated_lines[:start_idx_5] + new_block + updated_lines[end_idx:])

        # ===================BLOCK 8======================
        new_block = []
        new_block.append(f"{header_8}\n")
        new_block.append(f"{self.source_components:10d}\n")

        if self.source_components > 0:

            new_block.append(" ".join(f"{x:10.5f}" for x in self.source_lifetime) + "\n")
            new_block.append(" ".join(f"{x:10.3f}" for x in self.source_sigma) + "\n")
            new_block.append(" ".join(f"{x:10.5f}" for x in self.source_intensity) + "\n")
            new_block.append(f"{self.source_total:10.5f}\n")

            # Replace old block
            # We know structure length = 6 lines total
            end_idx = start_idx_8 + 6
        else:
            end_idx = start_idx_8 + 6 # still want to slice all 6 lines!

        updated_lines = (updated_lines[:start_idx_8] + new_block + updated_lines[end_idx:])

        with open(file_out, "w") as f:
            f.writelines(updated_lines)

    def run(self, exe_path, out_file=None, timeout=None):

        if self.path is None:
            raise ValueError("PFCFile.path is None.")

        exe_path = Path(exe_path)
        pfc_path = Path(self.path)

        if not exe_path.exists():
            raise FileNotFoundError(f"Executable not found: {exe_path}")

        if not pfc_path.exists():
            raise FileNotFoundError(f"PFC not found: {exe_path}")

        # Default output name in SAME directory as PFC
        if out_file is None:
            out_file = pfc_path.with_suffix(".out")

        out_file = Path(out_file)

        proc = subprocess.Popen([str(exe_path)], cwd=str(pfc_path.parent), stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,)

        # Only send filenames (not full paths)
        stdin_text = f"{pfc_path.name}\n{out_file.name}\n"

        stdout, stderr = proc.communicate(stdin_text, timeout=timeout)

        if proc.returncode != 0:
            raise RuntimeError(f"res19 failed (code {proc.returncode})\nSTDOUT:\n{stdout}\nSTDERR:\n{stderr}")

        produced = pfc_path.parent / out_file.name

        if not produced.exists():
            raise RuntimeError(f"Output file not created: {produced}")

        return produced

<<<<<<< Updated upstream
=======

class RFCOutFile:
    """Provides an OOP approach to the parameters within an .out file from a ResolutionFit (.rfc) run in PALSfit3
    A file can be read in with RFCOutFile.read()
    The attributes are:
        excursions
        job_time
        comment
        in_file
        dataset
        time_scale
        area_range
        fit_range
        init_fwhm
        init_intensity
        init_shifts
        converged
        iterations
        chi_square
        dof
        reduced_chi_square
        reduced_chi_square_std
        res_fwhm
        res_fwhm_std
        res_intensity
        res_intensity_std
        res_shift
        res_shift_std
        lt_lifetime
        lt_lifetime_std
        lt_intensity
        lt_intensity_std
        background
        background_std
        time_zero
        time_zero_std
        total_area_fit
        total_area_table
    """
    def __init__(self, path=None):
        self.path = Path(path) if path else None
        self._raw_lines = []

        # global fit stats
        self.converged = None
        self.iterations = None
        self.chi_square = None
        self.dof = None
        self.reduced_chi_square = None
        self.reduced_chi_square_std = None

        # Resolution final
        self.res_fwhm = None
        self.res_fwhm_std = None
        self.res_intensity = None
        self.res_intensity_std = None
        self.res_shift = None
        self.res_shift_std = None

        # Lifetime final
        self.lt_lifetime = None
        self.lt_lifetime_std = None
        self.lt_intensity = None
        self.lt_intensity_std = None

        # Background + time-zero
        self.background = None
        self.background_std = None
        self.time_zero = None
        self.time_zero_std = None

        # Total area
        self.total_area_fit = None
        self.total_area_table = None

    @classmethod
    def read(cls, file_in):
        obj = cls(file_in)
        file_in = Path(file_in)

        with open(file_in, "r") as f:
            obj._raw_lines = f.readlines()

        obj.parse_settings()
        return obj

    # -------- helpers --------
    _float_re = re.compile(r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?")

    @classmethod
    def _floats_in_line(cls, s):
        return [float(x) for x in cls._float_re.findall(s)]

    @classmethod
    def _stddevs_in_line_allow_fixed(cls, s):
        """
        Returns list containing floats and None for 'FIXED'.
        Assumes tokens after ':' are either 'FIXED' or numbers.
        """
        if ":" in s:
            s = s.split(":", 1)[1]
        toks = s.strip().split()
        out = []
        for t in toks:
            if t.upper().startswith("FIXED"):
                out.append(None)
            else:
                m = cls._float_re.fullmatch(t)
                if m:
                    out.append(float(t))
        return out

    # -------- parser --------
    def parse_settings(self):
        lines = self._raw_lines

        in_final = False
        in_res = False
        in_lt = False
        in_bg = False
        in_t0 = False

        for raw in lines:
            line = raw.rstrip("\n")

            # Enter final results section
            if "F I N A L  R E S U L T S" in line:
                in_final = True
                continue

            if not in_final:
                # still can parse some global lines if desired, but easiest is parse after FINAL
                continue

            # ---- global stats ----
            if "CONVERGENCE" in line and "ITERATIONS" in line:
                # Extract iteration count
                vals = self._floats_in_line(line)
                if vals:
                    self.iterations = int(vals[0])
                # Determine convergence status
                if "NOT OBTAINED" in line:
                    self.converged = False
                elif "OBTAINED" in line:
                    self.converged = True
                continue


            if "CHI-SQUARE" in line and "DEGREES OF FREEDOM" in line:
                vals = self._floats_in_line(line)
                if len(vals) >= 2:
                    self.chi_square = float(vals[0])
                    self.dof = int(vals[1])
                continue

            if "REDUCED CHI-SQUARE" in line:
                vals = self._floats_in_line(line)
                if len(vals) >= 2:
                    self.reduced_chi_square = float(vals[0])
                    self.reduced_chi_square_std = float(vals[1])
                continue

            # ---- section toggles ----
            if line.strip().startswith("RESOLUTION FUNCTION:"):
                in_res, in_lt, in_bg, in_t0 = True, False, False, False
                continue

            if line.strip().startswith("LIFETIME COMPONENTS:"):
                in_res, in_lt, in_bg, in_t0 = False, True, False, False
                continue

            if line.strip().startswith("BACKGROUND:"):
                in_res, in_lt, in_bg, in_t0 = False, False, True, False
                continue

            if line.strip().startswith("TIME-ZERO"):
                in_res, in_lt, in_bg, in_t0 = False, False, False, True
                # Note: the time-zero value is on the same line after ':'
                vals = self._floats_in_line(line)
                if vals:
                    self.time_zero = float(vals[0])
                continue

            # ---- parse inside RESOLUTION ----
            if in_res:
                if "FWHM (NS)" in line and ":" in line:
                    self.res_fwhm = self._floats_in_line(line)
                    continue
                if "INTENSITIES (%)" in line and ":" in line:
                    self.res_intensity = self._floats_in_line(line)
                    continue
                if "SHIFTS (NS)" in line and ":" in line:
                    self.res_shift = self._floats_in_line(line)
                    continue
                if "STD DEVIATIONS" in line and ":" in line:
                    # Which std dev line is it? Use what was last set but simplest:
                    # If res_fwhm exists and res_fwhm_std not yet set -> assign there, else if shift std not yet set -> assign there.
                    stds = self._stddevs_in_line_allow_fixed(line)
                    if self.res_fwhm is not None and self.res_fwhm_std is None:
                        self.res_fwhm_std = stds
                    elif self.res_shift is not None and self.res_shift_std is None:
                        self.res_shift_std = stds
                    continue

            # ---- parse inside LIFETIME ----
            if in_lt:
                if "LIFETIMES (NS)" in line and ":" in line:
                    self.lt_lifetime = self._floats_in_line(line)
                    continue
                if "INTENSITIES (%)" in line and ":" in line:
                    self.lt_intensity = self._floats_in_line(line)
                    continue
                if "STD DEVIATIONS" in line and ":" in line:
                    stds = self._stddevs_in_line_allow_fixed(line)
                    if self.lt_lifetime is not None and self.lt_lifetime_std is None:
                        self.lt_lifetime_std = stds
                    elif self.lt_intensity is not None and self.lt_intensity_std is None:
                        # intensity std devs are always numeric in your snippet
                        self.lt_intensity_std = [x for x in stds if x is not None]  # type: ignore
                    continue

            # ---- parse BACKGROUND ----
            if in_bg:
                if "COUNTS/CHANNEL" in line and ":" in line:
                    vals = self._floats_in_line(line)
                    if vals:
                        self.background = float(vals[0])
                    continue
                if "STD DEVIATION" in line and ":" in line:
                    vals = self._floats_in_line(line)
                    if vals:
                        self.background_std = float(vals[0])
                    continue

            # ---- parse TIME-ZERO std dev ----
            if in_t0:
                if "STD DEVIATIONS" in line and ":" in line:
                    vals = self._floats_in_line(line)
                    if vals:
                        self.time_zero_std = float(vals[0])
                    continue

            # ---- total area ----
            if line.strip().startswith("TOTAL AREA"):
                # "TOTAL AREA   FROM FIT         : 8.82842E+06     FROM TABLE : 7.20363E+06"
                vals = self._floats_in_line(line)
                if len(vals) >= 2:
                    self.total_area_fit = float(vals[0])
                    self.total_area_table = float(vals[1])
                continue


class PFCOutFile:
    def __init__(self, path=None):
        self.path = Path(path) if path else None
        self._raw_lines = []

        self.time_scale_ns_per_channel = None
        self.area_range_start_ch = None
        self.area_range_end_ch = None
        self.fit_range_start_ch = None
        self.fit_range_end_ch = None
        self.res_fwhm = None
        self.res_intensity = None
        self.res_shift = None

        # Before source correction
        self.no_corr_converged = None
        self.no_corr_iterations = None
        self.no_corr_chi_square = None
        self.no_corr_dof = None
        self.no_corr_reduced_chi_square = None
        self.no_corr_reduced_chi_square_std = None

        self.no_corr_lifetime = None
        self.no_corr_lifetime_std = None
        self.no_corr_intensity = None
        self.no_corr_intensity_std = None

        self.no_corr_background = None
        self.no_corr_background_std = None
        self.no_corr_time_zero = None
        self.no_corr_time_zero_std = None

        self.no_corr_total_area_fit = None
        self.no_corr_total_area_table = None

        # Source Correction
        self.source_lifetime = None
        self.source_intensity = None
        self.source_total = None

        # After source correction / Final results
        self.converged = None
        self.iterations = None
        self.chi_square = None
        self.dof = None
        self.reduced_chi_square = None
        self.reduced_chi_square_std = None
        self.lifetime = None
        self.lifetime_std = None
        self.sigma = None
        self.sigma_std = None
        self.intensity = None
        self.intensity_std = None
        self.mean_lifetime = None
        self.mean_lifetime_std = None
        self.background = None
        self.background_std = None
        self.time_zero = None
        self.time_zero_std = None
        self.total_area_fit = None
        self.total_area_table = None

    @classmethod
    def read(cls, file_in):
        obj = cls(file_in)
        file_in = Path(file_in)

        with open(file_in, "r") as f:
            obj._raw_lines = f.readlines()

        obj.parse_settings()
        return obj
    
    # -------- helpers --------
    _float_re = re.compile(r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?")

    @classmethod
    def _floats_in_line(cls, s):
        return [float(x) for x in cls._float_re.findall(s)]

    @classmethod
    def _stddevs_in_line_allow_fixed(cls, s):
        """
        Returns list containing floats and None for 'FIXED'.
        Assumes tokens after ':' are either 'FIXED' or numbers.
        """
        if ":" in s:
            s = s.split(":", 1)[1]
        toks = s.strip().split()
        out = []
        for t in toks:
            if t.upper().startswith("FIXED"):
                out.append(None)
            else:
                m = cls._float_re.fullmatch(t)
                if m:
                    out.append(float(t))
        return out
    
    # -------- parser --------
    def parse_settings(self):
        lines = self._raw_lines

        # ---- state flags ----
        in_initial = False
        in_no_corr_results = False
        in_source_corr = False
        in_final = False

        # ---- small helper ----
        def has(s, key):  # safe contains
            return key in s

        for raw in lines:
            line = raw.rstrip("\n")

            # -------------------------
            # Global / section switches
            # -------------------------
            if "----------------- I N I T I A L   P A R A M E T E R S ------------------" in line:
                in_initial = True
                in_no_corr_results = False
                in_source_corr = False
                in_final = False
                continue

            if "----- R E S U L T S  B E F O R E  S O U R C E  C O R R E C T I O N -----" in line:
                in_initial = False
                in_no_corr_results = True
                in_source_corr = False
                in_final = False
                continue

            if "------------------- S O U R C E  C O R R E C T I O N -------------------" in line:
                in_initial = False
                in_no_corr_results = False
                in_source_corr = True
                in_final = False
                continue

            # explicit "no source correction" banner (means: there is no before-corr results block)
            if "N O  S O U R C E  C O R R E C T I O N" in line:
                in_initial = False
                in_no_corr_results = False
                in_source_corr = False
                # final still coming; keep in_final False until we hit FINAL RESULTS banner
                continue

            if "####################### F I N A L  R E S U L T S #######################" in line:
                in_initial = False
                in_no_corr_results = False
                in_source_corr = False
                in_final = True
                continue

            # When we hit the end banner, stop parsing final
            if "######################### P O S I T R O N F I T ########################" in line:
                in_final = False
                continue

            # -------------------------
            # Parse "header-ish" values
            # -------------------------
            if has(line, "TIME SCALE") and ":" in line:
                vals = self._floats_in_line(line)
                if vals:
                    self.time_scale_ns_per_channel = float(vals[0])
                continue

            if has(line, "AREA RANGE") and "STARTS IN CH" in line and "ENDS IN CH" in line:
                vals = self._floats_in_line(line)
                if len(vals) >= 2:
                    self.area_range_start_ch = int(vals[0])
                    self.area_range_end_ch = int(vals[1])
                continue

            if has(line, "FIT RANGE") and "STARTS IN CH" in line and "ENDS IN CH" in line:
                vals = self._floats_in_line(line)
                if len(vals) >= 2:
                    self.fit_range_start_ch = int(vals[0])
                    self.fit_range_end_ch = int(vals[1])
                continue

            if has(line, "RESOLUTION") and has(line, "FWHM (NS)") and ":" in line:
                vals = self._floats_in_line(line)
                if vals:
                    self.res_fwhm = [float(x) for x in vals]  # can be 1 or 2 values
                continue

            if has(line, "FUNCTION") and has(line, "INTENSITIES") and ":" in line:
                vals = self._floats_in_line(line)
                if vals:
                    self.res_intensity = [float(x) for x in vals]
                continue

            if has(line, "SHIFTS (NS)") and ":" in line:
                vals = self._floats_in_line(line)
                if vals:
                    self.res_shift = [float(x) for x in vals]
                continue

            # -------------------------
            # Initial parameters block
            # -------------------------
            if in_initial:
                if line.strip().startswith("TIME-ZERO"):
                    vals = self._floats_in_line(line)
                    if vals:
                        self.init_time_zero = float(vals[0])
                    continue

                if line.strip().startswith("LIFETIMES (NS)"):
                    vals = self._floats_in_line(line)
                    if vals:
                        self.init_lifetime = [float(x) for x in vals]
                    continue

                if line.strip().startswith("SIGMA (NS)"):
                    vals = self._floats_in_line(line)
                    if vals:
                        self.init_sigma = [float(x) for x in vals]
                    continue

            # -----------------------------------------
            # Results BEFORE source correction (optional)
            # -----------------------------------------
            if in_no_corr_results:
                if "CONVERGENCE" in line and "ITERATIONS" in line:
                    vals = self._floats_in_line(line)
                    if vals:
                        self.no_corr_iterations = int(vals[0])
                    if "NOT OBTAINED" in line:
                        self.no_corr_converged = False
                    elif "OBTAINED" in line:
                        self.no_corr_converged = True
                    continue

                if "CHI-SQUARE" in line and "DEGREES OF FREEDOM" in line:
                    vals = self._floats_in_line(line)
                    if len(vals) >= 2:
                        self.no_corr_chi_square = float(vals[0])
                        self.no_corr_dof = int(vals[1])
                    continue

                # (some files may include reduced chi-square before-corr; your first example doesn't)
                if "REDUCED CHI-SQUARE" in line:
                    vals = self._floats_in_line(line)
                    if len(vals) >= 2:
                        self.no_corr_reduced_chi_square = float(vals[0])
                        self.no_corr_reduced_chi_square_std = float(vals[1])
                    continue

                if "LIFETIMES (NS)" in line and ":" in line:
                    vals = self._floats_in_line(line)
                    if vals:
                        self.no_corr_lifetime = [float(x) for x in vals]
                    continue

                if "SIGMA (NS)" in line and ":" in line:
                    vals = self._floats_in_line(line)
                    if vals:
                        self.no_corr_sigma = [float(x) for x in vals]
                    continue

                if "INTENSITIES (%)" in line and ":" in line:
                    vals = self._floats_in_line(line)
                    if vals:
                        self.no_corr_intensity = [float(x) for x in vals]
                    continue

                if line.strip().startswith("BACKGROUND"):
                    vals = self._floats_in_line(line)
                    if vals:
                        self.no_corr_background = float(vals[0])
                    continue

                if line.strip().startswith("TIME-ZERO"):
                    vals = self._floats_in_line(line)
                    if vals:
                        self.no_corr_time_zero = float(vals[0])
                    continue

                if line.strip().startswith("TOTAL AREA"):
                    # grabs both FROM FIT and FROM TABLE
                    vals = self._floats_in_line(line)
                    if len(vals) >= 2:
                        self.no_corr_total_area_fit = float(vals[0])
                        self.no_corr_total_area_table = float(vals[1])
                    continue

            # -------------------------
            # Source correction (optional)
            # -------------------------
            if in_source_corr:
                if "LIFETIMES (NS)" in line and ":" in line:
                    vals = self._floats_in_line(line)
                    if vals:
                        self.source_lifetime = [float(x) for x in vals]
                    continue

                if "INTENSITIES (%)" in line and ":" in line:
                    vals = self._floats_in_line(line)
                    if vals:
                        self.source_intensity = [float(x) for x in vals]
                    continue

                if "TOTAL (%)" in line and ":" in line:
                    vals = self._floats_in_line(line)
                    if vals:
                        self.source_total = float(vals[0])
                    continue

            # -------------------------
            # Final results (always)
            # -------------------------
            if in_final:
                if "CONVERGENCE" in line and "ITERATIONS" in line:
                    vals = self._floats_in_line(line)
                    if vals:
                        self.iterations = int(vals[0])
                    if "NOT OBTAINED" in line:
                        self.converged = False
                    elif "OBTAINED" in line:
                        self.converged = True
                    continue

                if "CHI-SQUARE" in line and "DEGREES OF FREEDOM" in line:
                    vals = self._floats_in_line(line)
                    if len(vals) >= 2:
                        self.chi_square = float(vals[0])
                        self.dof = int(vals[1])
                    continue

                if "REDUCED CHI-SQUARE" in line:
                    vals = self._floats_in_line(line)
                    if len(vals) >= 2:
                        self.reduced_chi_square = float(vals[0])
                        self.reduced_chi_square_std = float(vals[1])
                    continue

                # lifetimes and lifetime stddevs
                if "LIFETIMES (NS)" in line and ":" in line:
                    vals = self._floats_in_line(line)
                    if vals:
                        self.lifetime = [float(x) for x in vals]
                    continue

                if line.strip().startswith("STD DEVIATIONS") and (self.lifetime is not None) and (self.lifetime_std is None):
                    # first "STD DEVIATIONS" after lifetimes can include FIXED
                    vals = self._stddevs_in_line_allow_fixed(line)
                    if vals:
                        self.lifetime_std = vals
                    continue

                # sigma and sigma stddevs (can be *****)
                if "SIGMA (NS)" in line and ":" in line:
                    # if sigma line contains stars, floats_in_line may be empty; treat as None
                    vals = self._floats_in_line(line)
                    self.sigma = [float(x) for x in vals] if vals else None
                    continue

                if line.strip().startswith("STD DEVIATIONS") and ("sigma" in self.__dict__) and (self.sigma_std is None):
                    vals = self._stddevs_in_line_allow_fixed(line)
                    self.sigma_std = vals if vals else None
                    continue

                # intensities and intensity stddevs
                if "INTENSITIES (%)" in line and ":" in line:
                    vals = self._floats_in_line(line)
                    if vals:
                        self.intensity = [float(x) for x in vals]
                    continue

                if line.strip().startswith("STD DEVIATIONS") and (self.intensity is not None) and (self.intensity_std is None):
                    vals = self._stddevs_in_line_allow_fixed(line)
                    if vals:
                        self.intensity_std = vals
                    continue

                # mean lifetime + std
                if "MEAN LIFETIME" in line and ":" in line:
                    vals = self._floats_in_line(line)
                    if vals:
                        self.mean_lifetime = float(vals[0])
                    continue

                if line.strip().startswith("STD DEVIATION") and ("mean_lifetime" in self.__dict__) and (self.mean_lifetime_std is None):
                    vals = self._floats_in_line(line)
                    if vals:
                        self.mean_lifetime_std = float(vals[0])
                    continue

                # background + std
                if line.strip().startswith("BACKGROUND"):
                    vals = self._floats_in_line(line)
                    if vals:
                        self.background = float(vals[0])
                    continue

                if line.strip().startswith("STD DEVIATIONS") and (self.background is not None) and (self.background_std is None):
                    vals = self._floats_in_line(line)
                    if vals:
                        self.background_std = float(vals[0])
                    continue

                # time-zero + std
                if line.strip().startswith("TIME-ZERO"):
                    vals = self._floats_in_line(line)
                    if vals:
                        self.time_zero = float(vals[0])
                    continue

                if line.strip().startswith("STD DEVIATIONS") and (self.time_zero is not None) and (self.time_zero_std is None):
                    vals = self._floats_in_line(line)
                    if vals:
                        self.time_zero_std = float(vals[0])
                    continue

                # total area
                if line.strip().startswith("TOTAL AREA"):
                    vals = self._floats_in_line(line)
                    if len(vals) >= 2:
                        self.total_area_fit = float(vals[0])
                        self.total_area_table = float(vals[1])
                    continue
>>>>>>> Stashed changes
