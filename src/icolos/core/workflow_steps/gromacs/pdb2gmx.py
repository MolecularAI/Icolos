from icolos.utils.enums.program_parameters import (
    GromacsEnum,
)
from icolos.utils.enums.step_enums import StepGromacsEnum
from pydantic import BaseModel
from icolos.core.workflow_steps.gromacs.base import StepGromacsBase
from icolos.utils.execute_external.gromacs import GromacsExecutor
from icolos.utils.execute_external.execute import Executor
from icolos.utils.execute_external.schrodinger import SchrodingerExecutor
from icolos.core.workflow_steps.step import _LE
import os
import re
from typing import AnyStr, List
from string import ascii_uppercase

_SGE = StepGromacsEnum()
_GE = GromacsEnum()


class StepGMXPdb2gmx(StepGromacsBase, BaseModel):
    _shell_executor: Executor = None
    _antechamber_executor: Executor = None
    _acpype_executor: Executor = None
    _schrodinger_executor: SchrodingerExecutor = None

    def __init__(self, **data):
        """
        Executes system parametrisation for gromacs MD setup
        Generates GAFF params for unknown components with Antechamber
        """
        super().__init__(**data)

        self._initialize_backend(executor=GromacsExecutor)
        self._check_backend_availability()
        self._shell_executor = Executor()
        self._antechamber_executor = Executor()

    def _modify_topol_file(self, tmp_dir, itp_files):
        # read in the complex topol file, add the new itp files after the forcefield #include statement
        with open(os.path.join(tmp_dir, _SGE.COMPLEX_TOP), "r") as f:
            lines = f.readlines()
        index = [idx for idx, s in enumerate(lines) if _SGE.FORCEFIELD_ITP in s][0]
        new_topol = lines[: index + 1]
        for file in itp_files:
            new_topol.append(f'#include "{file}"\n')
        for line in lines[index + 1 :]:
            new_topol.append(line)
        for file in itp_files:
            stub = file.split(".")[0]
            new_topol.append(f"{stub}   1\n")
        with open(os.path.join(tmp_dir, _SGE.STD_TOPOL), "w") as f:
            f.writelines(new_topol)

        # remove all but the final topol file form the paths, makes file handling cleaner later
        top_files = [
            f for f in os.listdir(tmp_dir) if f.endswith("top") and f != _SGE.STD_TOPOL
        ]

        for f in top_files:
            os.remove(os.path.join(tmp_dir, f))

    def _add_posre_to_topol(self, tmp_dir, lig):
        """
        Add lines to topol file to invoke positional restraints for the parametrised ligands
        """
        stub = lig.split(".")[0]
        lig_itp = stub + ".itp"
        with open(os.path.join(tmp_dir, _SGE.STD_TOPOL), "r") as f:
            lines = f.readlines()
        index = [idx for idx, s in enumerate(lines) if lig_itp in s][0]
        new_topol = lines[: index + 1]
        new_topol.append(
            f"#ifdef POSRES_{stub.upper()}\n#include posre_{stub}.itp\n#endif\n"
        )
        for line in lines[index + 1 :]:
            new_topol.append(line)

        with open(os.path.join(tmp_dir, _SGE.STD_TOPOL), "w") as f:
            f.writelines(new_topol)

    def _split_protein_ligand_complex(self, tmp_dir):
        # split the file into protein and an arbitrary number of ligands and cofactors
        # Handle for multiple cofactors of the same type
        struct_file = [
            file for file in os.listdir(tmp_dir) if file.endswith(_SGE.FIELD_KEY_PDB)
        ][0]
        with open(os.path.join(tmp_dir, struct_file), "r") as f:
            data = f.readlines()
        # handles arbitrary number of ligands, cofactors, etc
        ligand_lines = {}
        protein_lines = []

        for line in data:
            parts = line.upper().split()

            # filter header lines etc
            if len(parts) > 4 and parts[0] in _GE.ATOMS:

                # catch the easy cases where there is a direct match to the parametrised components against internal dict
                if (
                    parts[3] in _GE.AMBER_PARAMETRISED_COMPONENTS
                    or parts[3] in _GE.IONS
                ):

                    protein_lines.append(line)

                # catch cases where ions have non-standard residue names e.g. NA3
                elif parts[3][:2] in _GE.IONS and re.findall(
                    re.compile(fr"{parts[3][:2]}[0-9]+"), line
                ):

                    pattern = fr"{parts[3][:2]}[0-9]+"
                    pattern = re.compile(pattern)

                    line = re.sub(pattern, parts[3][:2], line)
                    protein_lines.append(line)

                else:
                    # component is not parametrised, add to the ligands
                    if parts[4] in list(ascii_uppercase):
                        try:
                            ligand_lines[f"{parts[3]}:{parts[5]}"].append(line)
                        except KeyError:
                            # ligand key not created yet, identify by chain + res num to handle multiple identical components
                            ligand_lines[f"{parts[3]}:{parts[5]}"] = [line]
                    else:  # the 5th col index is the first coord col
                        try:
                            ligand_lines[f"{parts[3]}:{parts[4]}"].append(line)
                        except KeyError:
                            ligand_lines[f"{parts[3]}:{parts[4]}"] = [line]

        for key, value in ligand_lines.items():
            # write ligand components as separate pdb files
            with open(os.path.join(tmp_dir, f"{key}.pdb"), "w") as f:
                f.writelines(value)
        with open(os.path.join(tmp_dir, _SGE.PROTEIN_PDB), "w") as f:
            f.writelines(protein_lines)
        self._remove_temporary(os.path.join(tmp_dir, struct_file))
        return list(ligand_lines.keys())

    def _parametrisation_pipeline(self, tmp_dir, input_pdb) -> None:
        """
        :param tmp_dir: step's base directory
        :param input_pdb: file name for the ligand being parametrised
        """
        # main pipeline for producing GAFF parameters for a ligand
        stub = input_pdb.split(".")[0]
        output_file = stub + ".mol2"
        arguments_antechamber = [
            "-i",
            input_pdb,
            "-o",
            output_file,
            "-fi",
            "pdb",
            "-fo",
            "mol2",
            "-c",
            "gas",
        ]
        self._logger.log(f"Running antechamber on structure {input_pdb}", _LE.DEBUG)
        self._antechamber_executor.execute(
            command=_GE.ANTECHAMBER,
            arguments=arguments_antechamber,
            check=True,
            location=tmp_dir,
        )

        # Step 4: run the acpype script to generate the ligand topology file for GAFF
        charge_method = self.get_additional_setting(
            key=_SGE.CHARGE_METHOD, default="bcc"
        )
        self._logger.log(f"Running acpype on structure {input_pdb}", _LE.DEBUG)
        acpype_args = [
            "-di",
            output_file,
            "-c",
            charge_method,
        ]
        self._antechamber_executor.execute(
            command=_GE.ACPYPE_BINARY,
            arguments=acpype_args,
            location=tmp_dir,
            check=True,
        )
        # produce the ndx file for genrestr later
        index_file = stub + ".ndx"
        ndx_arguments = ["-f", input_pdb, "-o", index_file]

        self._backend_executor.execute(
            command=_GE.MAKE_NDX,
            arguments=ndx_arguments,
            location=tmp_dir,
            check=True,
            pipe_input='echo -e "0 & ! a H* \nq"',  # all system heavy atoms, excl hydrogens
        )
        # generate positional restraints for the ligand
        genrestr_args = [
            "-f",
            input_pdb,
            "-n",
            index_file,
            "-o",
            f"posre_{stub}.itp",
            "-fc",
            _SGE.FORCE_CONSTANTS,
        ]
        self._backend_executor.execute(
            command=_GE.GENRESTR,
            arguments=genrestr_args,
            location=tmp_dir,
            check=True,
            pipe_input="echo 3",
        )  # the ligand always be the last thing on the index file

        # we no longer need the ligand ndx file
        self._remove_temporary(os.path.join(tmp_dir, index_file))

    def _sort_components(self, lig_ids: List, components: List):
        """
        Ensure components go back into the concatenated pdb file in the same order as the original
        """
        new_components = []
        for idx in lig_ids:
            for component in components:
                if idx in component:
                    new_components.append(component)
        return new_components

    def _concatenate_structures(self, tmp_dir: str, lig_ids: List):
        """
        Extract newly parametrised components, concatenate everything into a single pdb file
        """

        components = []
        for root, _, files in os.walk(tmp_dir):
            for file in files:
                if file.endswith("_NEW.pdb"):
                    components.append(os.path.join(root, file))
        components = self._sort_components(lig_ids, components)
        self._logger.log(f"Found components: {components}", _LE.DEBUG)
        with open(os.path.join(tmp_dir, _SGE.PROTEIN_PDB), "r") as f:
            pdb_lines = f.readlines()

        for file in components:
            with open(file, "r") as f:

                pdb_lines.extend(f.readlines())

        pdb_lines = [
            l for l in pdb_lines if not any(s in l for s in ["TER", "ENDMDL", "REMARK"])
        ]
        pdb_lines.extend(["TER\n", "ENDMDL\n"])
        with open(os.path.join(tmp_dir, "Complex.pdb"), "w") as f:
            f.writelines(pdb_lines)

        # also deal with renaming the itp files here
        for root, _, files in os.walk(tmp_dir):
            for item in files:
                if (
                    item.endswith("GMX.itp")
                    and _SGE.PROTEIN_TOP not in item
                    and os.path.join(root, item) != os.path.join(tmp_dir, item)
                ):
                    os.rename(
                        os.path.join(root, item),
                        os.path.join(tmp_dir, item.split("_")[0]) + ".itp",
                    )
        # rename the protein top to complex
        os.rename(
            os.path.join(tmp_dir, _SGE.PROTEIN_TOP),
            os.path.join(tmp_dir, _SGE.COMPLEX_TOP),
        )

    def _extract_atomtype(self, tmp_dir: str, file: str) -> List[AnyStr]:
        """
        Pull the atomtype lines out of the topol file and return them as a list, write the sanitised itp file to directory
        """
        with open(os.path.join(tmp_dir, file), "r") as f:
            lines = f.readlines()
        start_index = None
        stop_index = None
        for idx, line in enumerate(lines):
            if _GE.ATOMTYPES in line:
                start_index = idx
            if _GE.MOLECULETYPES in line:
                stop_index = idx

        selection = lines[start_index:stop_index]
        # remove the offending lines from the topol
        remaining = lines[:start_index]
        remaining.extend(lines[stop_index:])
        self._remove_temporary(os.path.join(tmp_dir, file))
        with open(os.path.join(tmp_dir, file), "w") as f:
            f.writelines(remaining)
        return selection

    def _remove_duplicate_atomtypes(self, atomtypes: List):
        output = [atomtypes[0]]
        for line in atomtypes:
            if line not in output:
                output.append(line)
        return output

    def _modify_itp_files(self, tmp_dir):
        # cut the moleculetype directives out of all the individual itp files and add them to the top of the topol
        atomtype_lines = []
        # read the topol file, identify all the itp files it is #including
        with open(os.path.join(tmp_dir, _SGE.STD_TOPOL), "r") as f:
            topol_lines = [
                l.split()[-1].strip('"')
                for l in f.readlines()
                if ".itp" in l and "posre" not in l
            ]
        topol_lines = [l for l in topol_lines if l in os.listdir(tmp_dir)]
        for file in topol_lines:
            atomtype_lines.extend(self._extract_atomtype(tmp_dir, file))
        atomtype_lines = self._remove_duplicate_atomtypes(atomtype_lines)

        # write an 'atomtypes.itp' files to be included just below the forcefield, with all the atomtypes contained in the extra components
        with open(os.path.join(tmp_dir, "atomtypes.itp"), "w") as f:
            f.writelines(atomtype_lines)

        with open(os.path.join(tmp_dir, _SGE.STD_TOPOL), "r") as f:
            lines = f.readlines()
        self._remove_temporary(os.path.join(tmp_dir, _SGE.STD_TOPOL))
        index = [idx for idx, s in enumerate(lines) if _SGE.FORCEFIELD_ITP in s][0]
        new_topol = lines[: index + 1]

        new_topol.append('#include "atomtypes.itp"\n')
        new_topol.extend(lines[index + 1 :])
        with open(os.path.join(tmp_dir, _SGE.STD_TOPOL), "w") as f:
            f.writelines(new_topol)

    def _modify_water_molecules(self, tmp_dir: str):
        with open(os.path.join(tmp_dir, _SGE.COMPLEX_PDB), "r") as f:
            lines = f.readlines()

        solvent = []
        # pick out the water lines
        for line in lines:
            if any([x in line for x in _GE.SOLVENTS]):
                solvent.append(line)
        for line in solvent:
            lines.remove(line)
        lines.extend(solvent)
        for line in lines:
            if any([x in line for x in _GE.TERMINATIONS]):
                lines.remove(line)

        with open(os.path.join(tmp_dir, _SGE.COMPLEX_PDB), "w") as f:
            f.writelines(lines)

        if solvent:
            # modify the topol to put the solvent in last in the [ molecules ] directive
            with open(os.path.join(tmp_dir, _SGE.STD_TOPOL), "r") as f:
                lines = f.readlines()
            molecule_idx = lines.index(_GE.MOLECULES)
            for line in lines[molecule_idx:]:
                if any([x in line for x in _GE.SOLVENTS]):
                    out = lines.pop(lines.index(line))
                    lines.append(out)
            with open(os.path.join(tmp_dir, _SGE.STD_TOPOL), "w") as f:
                f.writelines(lines)

    def execute(self):
        """Takes in a ligand pdb file and generates the required topology, based on the backend specified in the config json.
        Currently supported AnteChamber

        Execution looks like this currently:
        (1) split the protein from the other components
        (2) generate the topology for the protein separately
        (3) identify components to be parametrised (cofactors, ligands etc)
        (4) run the parametrisation pipeline on each component in serial (reasonably fast exec time per ligand)
            (4a) store the resIDs of the ligands using the file handling system to be retrieved in a later step
        (5) modify the topology file to add the #include statements for the relevant itp files
        (6) convert the resulting concatenated pdb file to .gro with editconf
        (7) add the posres stuff to the topol file for each ligand for the subsequent equilibration steps
        (8) if more than one ligand, modify the itp files to ensure all moleculetype directives are specified first
        """

        tmp_dir = self._make_tmpdir()
        self._write_input_files(tmp_dir)  # dump generic data fields to the tmpdir
        lig_ids = self._split_protein_ligand_complex(tmp_dir)
        self._logger.log(
            f"Parameters will be generated for the following components: {str(lig_ids)}",
            _LE.DEBUG,
        )

        # Step 2: run pdb2gmx on the protein component only

        arguments_pdb2gmx = self._parse_arguments(
            flag_dict={
                "-f": os.path.join(tmp_dir, _SGE.PROTEIN_PDB),
                "-o": os.path.join(tmp_dir, _SGE.PROTEIN_PDB),
                "-p": _SGE.PROTEIN_TOP,
            }
        )
        self._backend_executor.execute(
            command=_GE.PDB2GMX, arguments=arguments_pdb2gmx, location=tmp_dir
        )

        for lig in lig_ids:
            input_file = lig + ".pdb"
            # generate the itp files for each component, named by their PDB identifier
            self._parametrisation_pipeline(tmp_dir, input_file)

        # concatenate the structures to produce Complex.pdb
        if lig_ids:
            self._concatenate_structures(tmp_dir, lig_ids)
            # step 6: Modify protein topol file for ligand
            itp_files = [
                f
                for f in os.listdir(tmp_dir)
                if f.endswith(".itp")
                and "posre" not in f
                and not any(
                    [x in f for x in _GE.PRIMARY_COMPONENTS]
                )  # avoid any duplicated itp file entries from components of the protein already handles by pdb2gmx (TODO: makes sure this works for DNA/RNA as well)
            ]
            # need to sort the itp files to match the ordering from the original pdb structure
            itp_files = self._sort_components(lig_ids, itp_files)
            self._modify_topol_file(tmp_dir, itp_files)

            # step 10: modify the topol file to add the ligand posre file if restraints are applied
            for lig in lig_ids:
                self._add_posre_to_topol(tmp_dir, lig)

            # if more than two ligands present, modify the ligand itp files so all the [atomtype] directives come before the [moleculetype] directives in the full topol
            if len(lig_ids) > 1:
                self._modify_itp_files(tmp_dir)

        else:
            # just convert the file names in place, no addition of ligands
            os.rename(
                os.path.join(tmp_dir, _SGE.PROTEIN_TOP),
                os.path.join(tmp_dir, _SGE.STD_TOPOL),
            )
            os.rename(
                os.path.join(tmp_dir, _SGE.PROTEIN_PDB),
                os.path.join(tmp_dir, _SGE.COMPLEX_PDB),
            )

            # step 7: run editconf to convert the combined pdb to a gro file

        # do final check to move crystallographic waters to the end of the pdb file, after
        # the ligand, to ensure continuous solvent group later
        self._modify_water_molecules(tmp_dir)
        # and adjust the topol file to put any solvent last

        editconf_arguments = ["-f", _SGE.COMPLEX_PDB, "-o", "structure.gro"]
        self._backend_executor.execute(
            command=_GE.EDITCONF,
            arguments=editconf_arguments,
            location=tmp_dir,
            check=True,
        )

        self._parse_output(tmp_dir)
        self._remove_temporary(tmp_dir)
