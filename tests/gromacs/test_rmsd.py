from icolos.core.containers.generic import GenericData
import unittest
import os
from icolos.core.workflow_steps.gromacs.rsmd import StepGMXrmsd
from icolos.utils.enums.step_enums import StepBaseEnum, StepGromacsEnum
from tests.tests_paths import (
    PATHS_EXAMPLEDATA,
    export_unit_test_env_vars,
    get_docked_ligands_as_conformers,
)
from icolos.utils.general.files_paths import attach_root_path

SGE = StepGromacsEnum()
SBE = StepBaseEnum


class Test_Trjcat(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/gromacs")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

        export_unit_test_env_vars()

    def setUp(self):
        with open(PATHS_EXAMPLEDATA.GROMACS_DMP_LIGAND_TRJ, "rb") as f:
            self.xtc = f.read()

        # load the docked pose as a compound
        self.comp = get_docked_ligands_as_conformers(
            PATHS_EXAMPLEDATA.GROMACS_DMP_LIGAND_SDF
        )
        print(self.comp)

    def test_gmx_rmsd(self):
        step_conf = {
            SBE.STEPID: "test_gmx_rmsd",
            SBE.STEP_TYPE: "gmx_rmsd",
            SBE.EXEC: {
                SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2021-fosscuda-2019a-PLUMED-2.7.1-Python-3.7.2"
            },
            SBE.SETTINGS: {},
        }

        step_rmsd = StepGMXrmsd(**step_conf)
        step_rmsd.data.generic.add_file(
            GenericData(file_name="structure.xtc", file_data=self.xtc, argument=True)
        )
        step_rmsd.data.compounds = self.comp

        step_rmsd.execute()
        out_path = os.path.join(self._test_dir, "rmsd.xvg")
        step_rmsd.write_generic_by_extension(self._test_dir, "xvg")
        stat_inf = os.stat(out_path)
        self.assertAlmostEqual(stat_inf.st_size, 3220, delta=100)
