import unittest
import os
from icolos.core.containers.generic import GenericData
from icolos.core.workflow_steps.openff.openff2gmx import StepOFF2gmx
from icolos.utils.enums.step_enums import StepBaseEnum, StepGromacsEnum, StepOpenFFEnum
from tests.tests_paths import (
    MAIN_CONFIG,
    PATHS_1UYD,
    PATHS_EXAMPLEDATA,
    create_test_dir,
    export_unit_test_env_vars,
)
from icolos.utils.general.files_paths import attach_root_path

_SBE = StepBaseEnum
_SOFE = StepOpenFFEnum()
_SGE = StepGromacsEnum()


class Test_PMXgenlib(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/openff/openff2gmx")

        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

        export_unit_test_env_vars()

    def setUp(self):
        with open(PATHS_1UYD.NATIVE_LIGAND_PDB, "r") as f:
            pdb_data = f.readlines()
        self.pdb_obj = GenericData(file_name="1UYD.pdb", file_data=pdb_data)

    def test_openff_ligand(self):
        step_conf = {
            _SBE.STEPID: "test_pdb2gmx",
            _SBE.STEP_TYPE: "openff2gmx",
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2021-fosscuda-2019a-PLUMED-2.7.1-Python-3.7.2"
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ADDITIONAL: {
                    _SOFE.UNIQUE_MOLS: [
                        "COc1c(OC)cc(c(Cl)c1OC)Cc(n2)n(CCCC)c(c23)ncnc3N"
                    ],
                    _SOFE.FORCEFIELD: os.path.join(
                        MAIN_CONFIG["OPENFF_FORCEFIELDS"], "openff-2.0.0.offxml"
                    ),
                },
            },
        }
        step_openff2gmx = StepOFF2gmx(**step_conf)
        step_openff2gmx.data.generic.add_file(self.pdb_obj)
        step_openff2gmx.execute()

        out_path = os.path.join(self._test_dir, "MOL.gro")
        step_openff2gmx.write_generic_by_extension(
            self._test_dir, _SGE.FIELD_KEY_STRUCTURE
        )
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 2400)
