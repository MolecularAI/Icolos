import shutil
from icolos.core.containers.gmx_state import GromacsState
from icolos.utils.enums.program_parameters import (
    GromacsEnum,
)
from icolos.utils.enums.step_enums import StepGromacsEnum, StepOpenFFEnum
from pydantic import BaseModel
from icolos.core.workflow_steps.gromacs.base import StepGromacsBase
from icolos.utils.execute_external.gromacs import GromacsExecutor
from icolos.utils.execute_external.execute import Executor
from icolos.utils.execute_external.schrodinger import SchrodingerExecutor
from icolos.core.workflow_steps.step import _LE
import os
import re
from string import ascii_uppercase
from rdkit import Chem
from openmm.app import PDBFile
from parmed.gromacs import GromacsTopologyFile
from openff.toolkit.topology import Molecule, Topology
from openff.toolkit.typing.engines.smirnoff import ForceField
import parmed as pmd

_SGE = StepGromacsEnum()
_GE = GromacsEnum()
_SOFE = StepOpenFFEnum()


class StepGMXPdb2gmx(StepGromacsBase, BaseModel):
    class Config:
        arbitrary_types_allowed = True

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
                    re.compile(rf"{parts[3][:2]}[0-9]+"), line
                ):

                    pattern = rf"{parts[3][:2]}[0-9]+"
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

    def _generate_gaff_params(
        self, tmp_dir: str, input_pdb: str, topol: GromacsState
    ) -> None:
        """
        :param tmp_dir: step's base directory
        :param input_pdb: file name for the ligand being parametrised

        Produces a gmx ITP file for the ligand, and updates topology
        """
        # main pipeline for producing GAFF parameters for a ligand
        charge_method = self.get_additional_setting(
            key=_SGE.CHARGE_METHOD, default="bcc"
        )
        conf = Chem.rdmolfiles.MolFromPDBFile(os.path.join(tmp_dir, input_pdb))
        formal_charge = Chem.rdmolops.GetFormalCharge(conf) if conf is not None else 0
        stub = input_pdb.split(".")[0]
        self._logger.log(
            f"Computed formal charge: {formal_charge} for structure {input_pdb}",
            _LE.DEBUG,
        )
        ligand_stub = input_pdb.split(".")[0]
        # Step 4: run the acpype script to generate the ligand topology file for GAFF
        self._logger.log(
            f"Running acpype on structure:{input_pdb}, charge method: {charge_method}",
            _LE.DEBUG,
        )
        acpype_args = [
            "-di",
            input_pdb,
            "-c",
            charge_method,
            "-o",
            "gmx",
            "-n",
            formal_charge,
        ]
        self._antechamber_executor.execute(
            command=_GE.ACPYPE_BINARY,
            arguments=acpype_args,
            location=tmp_dir,
            check=True,
        )

        acpype_dir = os.path.join(tmp_dir, f"{ligand_stub}.acpype")
        # TODO: refactor this
        lig_itp = [
            f
            for f in os.listdir(os.path.join(tmp_dir, acpype_dir))
            if f.endswith("GMX.itp")
        ][0]
        lig_pdb = [
            f
            for f in os.listdir(os.path.join(tmp_dir, acpype_dir))
            if f.endswith("pdb")
        ][0]

        topol.add_itp(os.path.join(tmp_dir, acpype_dir), [lig_itp], gen_posre=True)
        topol.add_molecule(lig_itp.split("_")[0], 1)
        topol.append_structure(
            os.path.join(tmp_dir, acpype_dir), lig_pdb, sanitize=True
        )

    def _generate_openff_params(
        self, tmp_dir: str, input_pdb: str, topol: GromacsState
    ):
        """
        Generate Amber-compatible Smirnoff params for a single component

        """
        # TODO: ensure we can handle parameter generation for multiple unique mol types
        # get the mols
        mols = [Molecule.from_smiles(self.get_additional_setting(_SOFE.UNIQUE_MOLS))]
        pdb_file = PDBFile(os.path.join(self.get_additional_setting("lig_pdb")))
        omm_topology = pdb_file.topology

        off_topology = Topology.from_openmm(omm_topology, unique_molecules=mols)

        forcefield = ForceField(self.get_additional_setting(_SOFE.FORCEFIELD))

        omm_system = forcefield.create_openmm_system(off_topology)

        parmed_struct = pmd.openmm.load_topology(
            omm_topology, system=omm_system, xyz=pdb_file.positions
        )
        parmed_struct.save(os.path.join(tmp_dir, "MOL.top"), overwrite=True)
        parmed_struct.save(os.path.join(tmp_dir, "MOL.pdb"), overwrite=True)
        gmx_top = GromacsTopologyFile().from_structure(parmed_struct)
        gmx_top.write(os.path.join(tmp_dir, "MOL.itp"), itp=True)
        topol.add_itp(path=tmp_dir, files=["MOL.itp"])
        topol.generate_posre(path=tmp_dir, itp_file="MOL.itp")
        topol.append_structure(path=tmp_dir, file="MOL.pdb", sanitize=True)

    def execute(self):
        """Takes in a ligand pdb file and generates the required topology, based on the backend specified in the config json.
        Currently supported AnteChamber

        Execution looks like this currently:
        (1) split the protein from the other components
        (2) generate the topology for the protein separately
        (4) run the parametrisation pipeline on each HET component in serial, generating either GAFF or SMIRNOFF params
        (5) combine the topologies

        """

        tmp_dir = self._make_tmpdir()
        topol = self.get_topol()
        self.write_input_files(tmp_dir)
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
        # instantiate a separate topology object that will do a beter job of writing out than parmed can

        topol.forcefield = self.settings.arguments.parameters["-ff"]
        topol.water = self.settings.arguments.parameters["-water"]
        topol.parse(tmp_dir, _SGE.PROTEIN_TOP)
        topol.set_structure(tmp_dir, _SGE.PROTEIN_PDB, sanitize=True)
        posre_files = [
            f for f in os.listdir(tmp_dir) if f.endswith(".itp") and "posre" in f
        ]
        topol.add_posre(tmp_dir, posre_files)

        param_method = self.get_additional_setting(_SGE.PARAM_METHOD, default=_SGE.GAFF)
        for lig in lig_ids:
            input_file = lig + ".pdb"
            if param_method == _SGE.GAFF:
                # generate the itp files for each component, named by their PDB identifier
                self._generate_gaff_params(tmp_dir, input_file, topol)
            elif param_method == _SGE.OPENFF:
                self._generate_openff_params(tmp_dir, input_file, topol)
            else:
                raise NotImplementedError

        os.remove(os.path.join(tmp_dir, _SGE.PROTEIN_TOP))
        os.remove(os.path.join(tmp_dir, _SGE.PROTEIN_PDB))
        # final writeout of parametrized system
        topol.write_structure(tmp_dir, _SGE.COMPLEX_PDB)
        # convert pdb to gro
        editconf_arguments = ["-f", _SGE.COMPLEX_PDB, "-o", _SGE.STD_STRUCTURE]
        self._backend_executor.execute(
            command=_GE.EDITCONF,
            arguments=editconf_arguments,
            location=tmp_dir,
            check=True,
        )
        topol.set_structure(tmp_dir)
        self._parse_output(tmp_dir)
        self._remove_temporary(tmp_dir)
