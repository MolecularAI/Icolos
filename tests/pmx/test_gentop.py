import unittest
import os
from icolos.core.composite_agents.workflow import WorkFlow
from icolos.core.containers.generic import GenericData
from icolos.core.containers.perturbation_map import Edge, PerturbationMap
from icolos.core.workflow_steps.pmx.gentop import StepPMXgentop
from icolos.utils.enums.step_enums import StepBaseEnum
from tests.tests_paths import PATHS_EXAMPLEDATA, export_unit_test_env_vars
from icolos.utils.general.files_paths import attach_root_path

_SBE = StepBaseEnum


class Test_PMXgentop(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/pmx/gentop")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

        export_unit_test_env_vars()

    def setUp(self):
        pass 

    def test_pmx_gentop(self):
        pass        