import unittest
import os

from icolos.core.containers.generic import GenericData, GenericContainer

from tests.tests_paths import PATHS_EXAMPLEDATA
from icolos.utils.general.files_paths import attach_root_path


class Test_Generic(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._GC = GenericContainer()

        cls._test_dir = attach_root_path("tests/junk/Generic")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    def setUp(self):
        gc = GenericContainer()
        with open(PATHS_EXAMPLEDATA.FEP_PLUS_DOCKING_PV, "rb") as f:
            data = f.read()
        gc.add_file(
            GenericData(file_name="test_file.txt", file_data=data, argument=True)
        )
        self.generic = gc

    def test_GenericHandling(self):
        self.assertEqual(len(self.generic.get_flattened_files()), 1)
        self.assertEqual(
            self.generic.get_file_by_name("test_file.txt").get_extension(), "txt"
        )
