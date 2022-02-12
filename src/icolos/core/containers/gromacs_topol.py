import os
from typing import AnyStr, Dict, List
from pydantic import BaseModel
from icolos.utils.enums.program_parameters import GromacsEnum

from icolos.utils.enums.step_enums import StepGromacsEnum

# acpype and pdb2gmx both do a job of handling posre creation etc

# this should be a singular object tacked on to the workflow that gets updated as the workflow progresses.

_SGE = StepGromacsEnum()
_GE = GromacsEnum()


class GromacsTopol(BaseModel):
    top_lines: List = []
    itps: Dict = {}
    posre: Dict = {}
    atomtypes: List = []
    chains: List = None
    forcefield: str = "amber03.ff"
    water: str = "tip3p"
    system: List = []
    molecules: Dict = {}
    # TODO: we could make this a generic data object?
    structure: List = []

    def __init__(self, **data) -> None:
        super().__init__(**data)

    def parse(self, file):
        """
        Populate the fields from parsing a topol file (normally from gmx pdb2gmx)
        Need system names and , molecules
        """
        # get the system name
        work_dir = os.path.dirname(file)
        with open(file, "r") as f:
            lines = f.readlines()
            # extract itp files (not including forcefield stuff)
        for line in lines:
            if line.startswith("#include") and ".itp" in line and ".ff" not in line:
                itp_file = line.split()[-1].replace('"', "")
                with open(os.path.join(work_dir, itp_file), "r") as f:
                    itp_lines = f.readlines()
                self.itps[itp_file] = itp_lines

        start = lines.index([l for l in lines if l.startswith(_GE.SYSTEM)][0])
        stop = lines.index([l for l in lines[start:] if l.startswith(_GE.MOLECULES)][0])
        for line in lines[start + 1 : stop]:
            if not line.startswith(";") and line.strip():
                self.system.append(line.strip())
        # excract molecules and counts
        for line in lines[stop + 1 :]:
            if not line.startswith(";") and line.strip():
                parts = line.split()
                self.molecules[parts[0].strip()] = parts[1].strip()

    def write_topol(self, base_dir: str, file: str = "topol.top"):
        """
        Write a gromacs topology file, including its dependent itp files, to a dir
        """
        if not self.top_lines:
            self.top_lines = self.generate_base_topol_file()
        with open(os.path.join(base_dir, file), "w") as f:
            for line in self.top_lines:
                f.write(line + "\n")
        for itp, lines in self.itps.items():
            if os.path.isfile(os.path.join(base_dir, itp)):
                os.remove(os.path.join(base_dir, itp))
            with open(os.path.join(base_dir, itp), "w") as f:
                for line in lines:
                    f.write(line)
        for itp, lines in self.posre.items():
            if os.path.isfile(os.path.join(base_dir, itp)):
                os.remove(os.path.join(base_dir, itp))
            with open(os.path.join(base_dir, itp), "w") as f:
                for line in lines:
                    f.write(line)

    def generate_base_topol_file(self) -> List[str]:
        """
        Generates the main topology file
        """
        top_lines = []
        top_lines.append(f'#include "{self.forcefield}.ff/forcefield.itp"')
        # find any new atomtypes in the parametrised components, slot these in first
        new_atomtypes = self.collect_atomtypes()
        top_lines += new_atomtypes

        # now include the itp files
        for file in self.itps.keys():
            top_lines.append(f'#include "{file}"')
            stub = file.split("_")[0]
            if "Protein" not in file:
                # these are handled independently in the itp files
                top_lines.append(f'#ifdef POSRES\n#include "posre_{stub}.itp"\n#endif')

        # add water model
        top_lines.append(f'#include "{self.forcefield}.ff/{self.water}.itp"')
        top_lines.append(_SGE.WATER_POSRE)

        # add ions itp
        top_lines.append(f'#include "{self.forcefield}.ff/ions.itp"')

        top_lines.extend(self._construct_block(_GE.SYSTEM, self.system))
        top_lines.extend(self._construct_block(_GE.MOLECULES, self.molecules))
        return top_lines

    def add_itp(self, path, files: List[str]) -> None:
        for file in files:
            with open(os.path.join(path, file), "r") as f:
                lines = f.readlines()
            self.itps[file] = lines

    def add_molecule(self, name: str, num: int = 1):
        self.molecules[name] = num

    def add_posre(self, path, files: List[str]) -> None:
        for file in files:
            with open(os.path.join(path, file), "r") as f:
                lines = f.readlines()

            self.posre[file] = lines

    def _construct_block(self, header, items) -> List:
        block = [header]
        if isinstance(items, dict):
            for key, value in items.items():
                block.append(f"{key}   {value}")
        else:
            for i in items:
                block.append(i)
        return block

    def collect_atomtypes(self) -> List:
        """
        Iterate over the itp files, strip any newly defined atomtypes, append to their own atomtypes.itp file, include this just after the forcefield include
        """
        atomtype_lines = []
        # read the topol file, identify all the itp files it is #including
        for itp, file_lines in self.itps.items():
            selection, remaining = self._extract_atomtype(file_lines)
            atomtype_lines.extend(selection)
            self.itps[itp] = remaining
        atomtype_lines = self._remove_duplicate_atomtypes(atomtype_lines)
        return atomtype_lines

    def _extract_atomtype(self, lines: List) -> List[AnyStr]:
        """
        Pull the atomtype lines out of the topol file and return them as a list, write the sanitised itp file to directory
        """

        start_index = None
        stop_index = None
        selection = []
        remaining = []
        for idx, line in enumerate(lines):
            if _GE.ATOMTYPES in line:
                start_index = idx
            if _GE.MOLECULETYPES in line:
                stop_index = idx
        if start_index is not None:
            selection = lines[start_index:stop_index]
            remaining = lines[:start_index]

        remaining.extend(lines[stop_index:])
        return selection, remaining

    def _remove_duplicate_atomtypes(self, atomtypes: List):
        output = [atomtypes[0]]
        for line in atomtypes:
            if line not in output:
                output.append(line)
        return output

    def set_topol(self, path: str, file: str = _SGE.STD_TOPOL):
        """
        When solvate or genion produce a modified topol, read this into the
        """
        with open(os.path.join(path, file), "r") as f:
            self.top_lines = f.readlines()

    def set_structure(self, path: str, file: str = _SGE.STD_STRUCTURE, sanitize=False):
        with open(os.path.join(path, file), "r") as f:
            lines = f.readlines()
        if sanitize:
            lines = [l for l in lines if any([l.startswith(idx) for idx in _GE.ATOMS])]
        self.structure = lines

    def write_structure(self, path: str, file: str = _SGE.STD_STRUCTURE):
        with open(os.path.join(path, file), "w") as f:
            for line in self.structure:
                f.write(line)

    def append_structure(self, path: str, file: str, sanitize=False) -> None:
        with open(os.path.join(path, file), "r") as f:
            lines = f.readlines()
        if sanitize:
            lines = [l for l in lines if any([l.startswith(idx) for idx in _GE.ATOMS])]
        self.structure.extend(lines)

    def __str__(self) -> str:
        return f"Gromacs Topology object: System: {self.system} | Molecules: {[m for m in self.molecules]} | FF: {self.forcefield} | itp files: {[f for f in self.itps.keys()]}"

    def __repr__(self) -> str:
        return self.__str__()
