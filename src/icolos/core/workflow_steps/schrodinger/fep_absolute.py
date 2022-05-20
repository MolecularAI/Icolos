import os
from pydantic import BaseModel
from icolos.core.workflow_steps.schrodinger.fep_base import StepFEPBase
from icolos.utils.enums.logging_enums import LoggingConfigEnum
from icolos.utils.enums.program_parameters import (
    FepPlusEnum,
    SchrodingerExecutablesEnum,
)
from icolos.utils.enums.write_out_enums import WriteOutEnum
from icolos.utils.execute_external.fep_plus import FepPlusExecutor
from rdkit import Chem

from icolos.utils.execute_external.schrodinger import SchrodingerExecutor

_WE = WriteOutEnum()
_SEE = SchrodingerExecutablesEnum()
_FE = FepPlusEnum()
_LE = LoggingConfigEnum()


class StepSchrodingerAbsoluteFEP(StepFEPBase, BaseModel):
    _schrodinger_executor: SchrodingerExecutor = None

    def __init__(self, **data):
        super().__init__(**data)
        self._initialize_backend(executor=FepPlusExecutor)
        self._check_backend_availability()

        self._schrodinger_executor = SchrodingerExecutor(
            prefix_execution=self.execution.prefix_execution,
            binary_location=self.execution.binary_location,
        )

    def _prepare_poseviewer(self, tmp_dir: str) -> None:

        # we must have a protein structure to construct the pv by combining with existing compounds
        self.data.generic.get_argument_by_extension("pdb", rtn_file_object=True).write(
            os.path.join(tmp_dir, "receptor.pdb"), join=False
        )

        # write conformers in self.data.compound to sdf
        with Chem.SDWriter(os.path.join(tmp_dir, "compounds.sdf")) as writer:
            for compound in self.get_compounds():
                for enum in compound.get_enumerations():
                    for conf in enum.get_conformers():
                        mol = conf.get_molecule()
                        mol_name = (
                            compound.get_name()
                            if compound.get_name()
                            else compound.get_index_string()
                        )
                        mol.SetProp(_WE.RDKIT_NAME, mol_name)
                        writer.write(mol)

        # run structcat with the ligand files
        structcat_args = [
            "-ipdb",
            "receptor.pdb",
            "-isd",
            "compounds.sdf",
            "-omae",
            "abfe_pv.maegz",
        ]
        self._schrodinger_executor.execute(
            command=_SEE.STRUCTCAT,
            arguments=structcat_args,
            check=True,
            location=tmp_dir,
        )
        # schrodinger return codes not always reliable
        assert os.path.isfile(os.path.join(tmp_dir, "abfe_pv.maegz"))

    def execute(self):
        """Execute Schrodinger's ABFE workflow"""

        tmp_dir = self._make_tmpdir()
        print(tmp_dir)

        # require poseviewer file, either direct from generic data or construct from receptor and ligands from previous dock
        if self.data.generic.get_files_by_extension("maegz"):
            # poseviewer has been provided
            self.data.generic.get_argument_by_extension("maegz").write(
                os.path.join(tmp_dir, "abfe_pv.maegz"), join=False
            )
            self._logger.log("Using provided poseviewer file", _LE.DEBUG)
        else:
            self._logger.log(
                "Preparing poseviewer from receptor and compounds", _LE.DEBUG
            )
            self._prepare_poseviewer(tmp_dir)

        abfe_args = ["abfe_pv.maegz"] + self._parse_arguments()
        self._backend_executor.execute(
            command=_FE.FEP_ABSOLUTE_EXECUTOR,
            arguments=abfe_args,
            location=tmp_dir,
            check=True,
        )
        self._parse_output(tmp_dir)
        self._remove_temporary(tmp_dir)
