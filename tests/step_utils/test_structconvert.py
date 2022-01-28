import os
import unittest
from icolos.core.step_utils.structconvert import StructConvert
from icolos.utils.enums.program_parameters import SchrodingerExecutablesEnum
from icolos.utils.enums.step_enums import StepBaseEnum
from icolos.utils.general.files_paths import attach_root_path, remove_folder
from tests.tests_paths import PATHS_EXAMPLEDATA

_SBE = StepBaseEnum
_SEE = SchrodingerExecutablesEnum()


class Test_Structconvert(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/structconvert")
        remove_folder(cls._test_dir)
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    def setUp(self):
        pass

    def test_sdf2pdb(self):
        executor = StructConvert(prefix_execution=_SEE.SCHRODINGER_MODULE)
        output_path = os.path.join(self._test_dir, "output_small_molecule.pdb")
        executor.sdf2pdb(
            sdf_file=PATHS_EXAMPLEDATA.SMALL_MOLECULE_SDF_PATH, pdb_file=output_path
        )

        stat_inf = os.stat(output_path)
        self.assertEqual(stat_inf.st_size, 2209)
