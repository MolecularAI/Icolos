import unittest
from icolos.core.workflow_steps.schrodinger.fep_absolute import (
    StepSchrodingerAbsoluteFEP,
)

from icolos.utils.enums.step_enums import StepBaseEnum
from icolos.core.containers.generic import GenericData
from tests.tests_paths import (
    MAIN_CONFIG,
    PATHS_1UYD,
    PATHS_EXAMPLEDATA,
)
import os
from icolos.utils.enums.step_enums import StepBaseEnum, StepGromacsEnum, StepPrepwizEnum
from tests.tests_paths import (
    PATHS_EXAMPLEDATA,
    get_ligands_as_compounds_with_conformers,
)
from icolos.utils.general.files_paths import attach_root_path
from icolos.core.workflow_steps.schrodinger.prepwizard import StepPrepwizard

_SGE = StepGromacsEnum()
_SBE = StepBaseEnum
_SPE = StepPrepwizEnum()


class Test_SchrodingerABFE(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/schrodinger_abfe")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    def setUp(self):
        with open(PATHS_1UYD.PDB_PATH, "r") as f:
            apo_1uyd = f.read()
        self.apo_1uyd = GenericData(file_name="receptor.pdb", file_data=apo_1uyd)
        self.compounds = get_ligands_as_compounds_with_conformers(PATHS_1UYD.LIGANDS)

    def test_abfe_construct_pv(self):
        step_conf = {
            _SBE.STEPID: "01_abfe",
            _SBE.STEP_TYPE: _SBE.STEP_FEP_ABSOLUTE,
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load schrodinger/2021-2-js-aws",
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                        "-HOST": "cpu-only",
                        "-SUBHOST": "fep-compute",
                        "-ppj": "1",
                        "-JOBNAME": "abfe_test_job",
                    }
                }
            },
        }

        abfe_step = StepSchrodingerAbsoluteFEP(**step_conf)
        abfe_step.data.generic.add_file(self.apo_1uyd)
        abfe_step.data.compounds = self.compounds
        # FIXME: this is a very  long execution, do not run as part of normal test suite!
        # abfe_step.execute()

        # out_file = abfe_step.data.generic.get_files_by_extension("pdb")[0].get_data()
        # out_path = os.path.join(self._test_dir, "test_out.pdb")
        # with open(out_path, "w") as f:
        #     f.write(out_file)
        # stat_inf = os.stat(out_path)
        # self.assertGreater(stat_inf.st_size, 155800)
