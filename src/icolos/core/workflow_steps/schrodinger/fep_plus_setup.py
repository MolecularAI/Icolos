from typing import List
from icolos.core.containers.generic import GenericData
from icolos.core.step_utils.structconvert import StructConvert
from icolos.core.workflow_steps.schrodinger.base import StepSchrodingerBase
from icolos.utils.enums.program_parameters import (
    FepPlusEnum,
    SchrodingerExecutablesEnum,
)
from icolos.utils.enums.step_enums import StepBaseEnum, StepFepPlusEnum, StepGlideEnum
from icolos.utils.execute_external.fep_plus import FepPlusExecutor
from rdkit.Chem import SDMolSupplier
from icolos.utils.execute_external.schrodinger import SchrodingerExecutor
from icolos.core.workflow_steps.step import _LE
import os
from pydantic import BaseModel
from rdkit.Chem import SDWriter

_SFE = StepFepPlusEnum()
_FE = FepPlusEnum()
_SEE = SchrodingerExecutablesEnum()
_SBE = StepBaseEnum
_SGE = StepGlideEnum()


class StepFepPlusSetup(StepSchrodingerBase, BaseModel):
    """
    Construct and analyse perturbation map for set of congeneric ligands
    Supports extracting structures from poseviewer or pdb files
    """

    _schrodinger_executor: SchrodingerExecutor = None
    _converter: StructConvert = None

    def __init__(self, **data):
        super().__init__(**data)

        self._initialize_backend(executor=FepPlusExecutor)
        self._check_backend_availability()

        self._schrodinger_executor = SchrodingerExecutor(
            prefix_execution=self.execution.prefix_execution,
            binary_location=self.execution.binary_location,
        )
        self._converter = StructConvert(
            prefix_execution=self.execution.prefix_execution,
            binary_location=self.execution.binary_location,
        )

    def _extract_receptor_from_pv(self, tmp_dir, input_file: str = _SFE.RECEPTOR_MAEGZ):
        # run split_structure.py to obtain the receptor_structure
        self._logger.log(f"Extracting receptor from structure.", _LE.INFO)
        self._schrodinger_executor.execute(
            command=_SEE.STRUCT_SPLIT,
            arguments=[
                "-m",
                "pdb",
                "-many_files",
                os.path.join(tmp_dir, input_file),
                f"{_SFE.STRUCT_SPLIT_BASE}.mae",
            ],
            check=True,
            location=tmp_dir,
        )

        # get rid of the original receptor structure now we have the new one
        os.remove(os.path.join(tmp_dir, _SFE.RECEPTOR_MAEGZ))

    def _write_receptor_from_pv(self, tmp_dir):
        # Handles writing the receptor structure to tmpdir, either from a poseviewer file, or a provided receptor
        # take the first poseviewer file it can find and split the stricure, take only the receptor
        for compound in self.get_compounds():
            for enumeration in compound.get_enumerations():
                for conformer in enumeration.get_conformers():
                    if (
                        _SGE.GLIDE_POSEVIEWER_FILE_KEY
                        in conformer.get_extra_data().keys()
                    ):
                        with open(
                            os.path.join(tmp_dir, _SFE.RECEPTOR_MAEGZ), "wb"
                        ) as f:
                            f.write(
                                conformer.get_extra_data()[
                                    _SGE.GLIDE_POSEVIEWER_FILE_KEY
                                ]
                            )
                            break
        if _SFE.RECEPTOR_MAEGZ in os.listdir(tmp_dir):
            self._logger.log(
                f"Writing poseviewer file to temporary directory.", _LE.INFO
            )
            self._extract_receptor_from_pv(tmp_dir)
        elif self.data.generic.get_files_by_extension("pdb"):
            # a pdb file was loaded to generic data, use this as the receptor structure
            self.data.generic.get_argument_by_extension(
                "pdb", rtn_file_object=True
            ).write(os.path.join(tmp_dir, "receptor.pdb"), join=False)

            self._logger.log(
                "Converting provided pdb receptor structure to mae", _LE.DEBUG
            )
            self._converter.convert(
                os.path.join(tmp_dir, "receptor.pdb"),
                os.path.join(tmp_dir, f"{_SFE.STRUCT_SPLIT_BASE}_receptor1.mae"),
            )
            os.remove(os.path.join(tmp_dir, "receptor.pdb"))

        else:
            self._logger.log(
                "No poseviewer file was found attached to any of the conformers, and no PDB receptor file was specified - this must be set in the docking step",
                _LE.ERROR,
            )
            raise FileNotFoundError

    def _check_xray_structure(self, compound_number):
        # check to see if an xray structure has been provided for that compound
        if _SFE.XRAY_STRUCTURES in self.settings.additional.keys():
            if isinstance(self.settings.additional[_SFE.XRAY_STRUCTURES], dict):
                if (
                    compound_number
                    in self.settings.additional[_SFE.XRAY_STRUCTURES].keys()
                ):
                    return True, _FE.DICT
            elif os.path.isdir(self.settings.additional[_SFE.XRAY_STRUCTURES]):
                if os.path.isfile(
                    os.path.join(
                        self.settings.additional[_SFE.XRAY_STRUCTURES],
                        f"{compound_number}.pdb",
                    )
                ):
                    return True, _FE.PATH
        return False, None

    def _rename_sdf(self, path, comp_num):
        with open(path, "r") as f:
            lines = f.readlines()[1:]
        new_lines = [f"{comp_num}:0:0\n"]
        for line in lines:
            new_lines.append(line)
        self._remove_temporary(path)
        with open(path, "w") as f:
            f.writelines(new_lines)

    def _extract_ligand_from_pdb(self, tmp_dir: str, comp_num: int, type: str):
        # if ligand poses have been provided from xray structures, extract just the ligand
        self._logger.log(
            f"Extracting ligand from provided Xray structure for compound {comp_num}",
            _LE.DEBUG,
        )
        if type == _FE.DICT:
            file_path = self.settings.additional[_SFE.XRAY_STRUCTURES[comp_num]]
        else:
            file_path = os.path.join(
                self.settings.additional[_SFE.XRAY_STRUCTURES], f"{comp_num}.pdb"
            )
        if not os.path.isfile(file_path):
            raise FileNotFoundError(
                "The provided path to the xray structure does not exist or is not accessible"
            )
        self._schrodinger_executor.execute(
            command=_SEE.STRUCT_SPLIT,
            arguments=["-m", "pdb", "-many_files", file_path, f"{_SFE.XRAY_SPLIT}.sdf"],
            check=True,
            location=tmp_dir,
        )
        # remove everything apart from the ligand sdf which is concatenated later
        lig_found = False
        for file in os.listdir(tmp_dir):
            idx = file.split("/")[-1]
            if idx.startswith(_SFE.XRAY_SPLIT):
                if "ligand" in idx:
                    # need to modify the name from the standard that Schrodinger provides
                    self._rename_sdf(os.path.join(tmp_dir, file), comp_num)
                    mols = SDMolSupplier(os.path.join(tmp_dir, file))

                    data = mols[0]
                    lig_found = True
                    self._remove_temporary(os.path.join(tmp_dir, file))
                else:
                    self._remove_temporary(os.path.join(tmp_dir, file))
        if lig_found:
            return data

    def _write_input_files(self, tmp_dir):
        # write receptor structure to tmpdir, either from poseviewer or provided pdb file
        self._write_receptor_from_pv(tmp_dir)

        # write out all conformers present in self.data.compounds to a single sdf file.
        writer = SDWriter(os.path.join(tmp_dir, "concatenated.sdf"))
        for compound in self.get_compounds():
            # If an xray pose is provided, use this
            flag, type = self._check_xray_structure(compound.get_compound_number())
            if flag is True:
                self._logger.log(
                    "Found Xray structure for the ligand - using this in preference to a docking pose",
                    _LE.DEBUG,
                )
                mol = self._extract_ligand_from_pdb(
                    tmp_dir, compound.get_compound_number(), type
                )
                writer.write(mol)
            else:
                # use the docked conformer
                for enumeration in compound.get_enumerations():
                    for conformer in enumeration.get_conformers():
                        mol = conformer.get_molecule()
                        writer.write(mol)

    def _parse_arguments(self, io_dict: dict) -> List[str]:
        arguments = []
        for key in self.settings.arguments.parameters.keys():
            arguments.append(key)
            arguments.append(str(self.settings.arguments.parameters[key]))
        for flag in self.settings.arguments.flags:
            arguments.append(str(flag))
        for key, value in io_dict.items():
            arguments.append(key)
            arguments.append(value)
        return arguments

    def _get_structcat_args(
        self, tmp_dir: str, out_file_type: str, outfile: str
    ) -> List[str]:
        arguments = [
            f"{_SEE.STRUCTCAT_I}mae",
            os.path.join(tmp_dir, f"{_SFE.STRUCT_SPLIT_BASE}_receptor1.mae"),
            f"{_SEE.STRUCTCAT_I}sd",
        ]

        for file in os.listdir(tmp_dir):
            if file.endswith("sdf"):
                arguments.append(os.path.join(tmp_dir, file))
        arguments.append(f"{_SEE.STRUCTCAT_O}{out_file_type}")
        arguments.append(os.path.join(tmp_dir, outfile))
        return arguments

    def _concatenate_pv_files(self, tmp_dir: str):
        # create a poseviewer-formatted file with receptor structure, then docked ligand poses
        arguments = self._get_structcat_args(
            tmp_dir=tmp_dir, out_file_type="mae", outfile=_SFE.STRUCTCAT_MAEGZ_OUTFILE
        )
        self._schrodinger_executor.execute(
            command=_SEE.STRUCTCAT, arguments=arguments, check=True
        )

    def _analyse_map(self, tmp_dir):
        """run fmp_stats program to analyse map - generate node similarities etc"""
        result = self._schrodinger_executor.execute(
            command=_SEE.FMP_STATS,
            arguments=["out.fmp", "-f"],
            check=True,
            location=tmp_dir,
        )
        log_lines = []
        for line in str(result.stdout).split("\n"):
            self._logger_blank.log(line, _LE.INFO)
            log_lines.append(line + "\n")

        self.data.generic.add_file(
            GenericData(file_name="fep_mapper.log", file_data=log_lines)
        )

    def _parse_output(self, tmp_dir: str):
        # needs to retrieve the edge and fmp files produced by the mapper step and attach to the generic dict
        files = [
            os.path.join(tmp_dir, f)
            for f in os.listdir(tmp_dir)
            if f.endswith(("fmp", "edge", "log"))
        ]

        for file in files:
            try:
                with open(file, "r") as f:
                    data = f.read()
            except UnicodeDecodeError:
                with open(file, "rb") as f:
                    data = f.read()
            self._add_data_to_generic(file, data)

    def execute(self):
        # run the job in a temporary directory
        tmp_dir = self._make_tmpdir()

        self._write_input_files(tmp_dir)
        self._concatenate_pv_files(tmp_dir)
        io_dict = {
            "": os.path.join(tmp_dir, _SFE.STRUCTCAT_MAEGZ_OUTFILE),
            "-o": _SFE.FEP_MAPPER_OUTPUT,
        }
        arguments = self._parse_arguments(io_dict=io_dict)
        self._apply_token_guard()  # need to implement for reliability
        self._logger.log("Optimising perturbation map", _LE.DEBUG)
        self._backend_executor.execute(
            command=_FE.FEP_MAPPER, arguments=arguments, check=True, location=tmp_dir
        )
        assert os.path.isfile(os.path.join(tmp_dir, "out.fmp"))
        self._logger.log(
            f"Successfully executed fep_mapper in directory {tmp_dir}.", _LE.DEBUG
        )

        self._logger.log("Analysing the perturbation map.", _LE.DEBUG)
        self._analyse_map(tmp_dir)
        self._parse_output(tmp_dir)
        self._remove_temporary(tmp_dir)
