import unittest
import os
from tests.tests_paths import export_unit_test_env_vars
from icolos.utils.general.files_paths import attach_root_path


class Test_PMXdoublebox(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/pmx/doublebox")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

        export_unit_test_env_vars()

    def setUp(self):
        pass

    def test_XYZ(self):
        pass
