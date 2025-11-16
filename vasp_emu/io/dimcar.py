"""Handles formatting and writing of OUTCAR"""
from ase.atoms import Atoms
from logging import Logger, FileHandler, Formatter
import numpy as np


class DimcarWriter(Logger):
    """
    A class to handle writing to the DIMCAR file.
    Inherits from Logger to use logging functionality.
    """
    def __init__(self, dir: str = ".", isif: int = 0) -> None:
        """
        Initialize the DimcarWriter at `dir/dimcar`.

        Args:
            dir (str): Directory where the OUTCAR file will be written.
        """
        super().__init__(name=dir + "/dimcar")
        # Prevent adding multiple handlers if the logger instance is somehow reused by name
        if not self.handlers:
            self.propagate = False
            self.setLevel("INFO")
            file_handler = FileHandler(filename=f"{dir}/DIMCAR", mode="w")
            file_handler.setFormatter(Formatter("%(message)s"))
            self.addHandler(file_handler)
            self.isif = isif

    def write_dimcar_header(self) -> None:
        """
        writes the header for dimcar
        """
        # Write the header for the file
        self.info(" Step         Force        Torque        Energy     Curvature         Angle")\

    def write_dimcar(self, data: dict) -> None:
        """
        Takes the dimcar data and writes the dimcar

        Args:
            data (data): should include step#, fmax, torque, 
                         curvature, and angles per step
        """
        
        # Write all the angles for each step
        for i in range(len(data["angles"])):
            self.info(f"{data['step']:>5}{data['fmax']:>14.5f}{data['torques'][i]:>14.5f}{data['energy']:>14.5f}{data['curvatures'][i]:>14.5f}{data['angles'][i]:>14.5f}")

