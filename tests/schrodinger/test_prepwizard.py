import unittest

from icolos.utils.enums.step_enums import StepBaseEnum
from icolos.core.workflow_steps.schrodinger.prepwizard import StepPrepwizard
from icolos.core.containers.generic import GenericData
from tests.tests_paths import (
    PATHS_1UYD,
    PATHS_EXAMPLEDATA,
)
import os
from icolos.utils.enums.step_enums import StepBaseEnum, StepGromacsEnum, StepPrepwizEnum
from tests.tests_paths import PATHS_EXAMPLEDATA
from icolos.utils.general.files_paths import attach_root_path
from icolos.core.workflow_steps.schrodinger.prepwizard import StepPrepwizard

_SGE = StepGromacsEnum()
_SBE = StepBaseEnum
_SPE = StepPrepwizEnum()


class Test_Prepwizard(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/prepwizard")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    def setUp(self):
        with open(PATHS_1UYD.PDB_PATH, "r") as f:
            apo_1uyd = f.read()
        self.apo_1uyd = GenericData(file_name="test_structure.pdb", file_data=apo_1uyd)
        with open(PATHS_EXAMPLEDATA.GROMACS_HOLO_STRUCTURE, "r") as f:
            self.holo_1uyd = f.read()

    def test_prepwizard(self):
        step_conf = {
            _SBE.STEPID: "01_ligprep",
            _SBE.STEP_TYPE: _SBE.STEP_PREPWIZARD,
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load schrodinger/2020-4",
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {_SBE.SETTINGS_ARGUMENTS_PARAMETERS: {}}
            },
        }

        prepwiz_step = StepPrepwizard(**step_conf)
        prepwiz_step.data.generic.add_file(self.apo_1uyd)
        prepwiz_step.execute()

        out_file = prepwiz_step.data.generic.get_files_by_extension("pdb")[0].get_data()
        out_path = os.path.join(self._test_dir, "test_out.pdb")
        with open(out_path, "w") as f:
            f.write(out_file)
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 155800)

    def test_remove_ligand(self):
        step_conf = {
            _SBE.STEPID: "test_rem",
            _SBE.STEP_TYPE: _SBE.STEP_PREPWIZARD,
            _SBE.EXEC: {_SBE.EXEC_PREFIXEXECUTION: "ml schrodinger"},
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ADDITIONAL: {_SPE.REMOVE_RES: ["DMP"]},
            },
        }

        step_removelig = StepPrepwizard(**step_conf)
        step_removelig.data.generic.add_file(
            GenericData(file_name="1BVG.pdb", file_data=self.holo_1uyd, argument=True)
        )

        step_removelig.execute()
        out_path = os.path.join(self._test_dir, "cox.pdb")
        step_removelig.write_generic_by_extension(
            self._test_dir,
            _SGE.PROTEIN_PDB,
        )

        out_file = step_removelig.data.generic.get_files_by_extension("pdb")[
            0
        ].get_data()
        with open(out_path, "w") as f:
            f.write(out_file)
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 126530)

    def test_auto_remove_ligand(self):
        step_conf = {
            _SBE.STEPID: "test_rem",
            _SBE.STEP_TYPE: _SBE.STEP_PREPWIZARD,
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load schrodinger/2021-2-js-aws"
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ADDITIONAL: {_SPE.REMOVE_RES: "ligands"},
            },
        }

        step_removelig = StepPrepwizard(**step_conf)
        step_removelig.data.generic.add_file(
            GenericData(
                file_name="1UYD_holo.pdb", file_data=self.holo_1uyd, argument=True
            )
        )

        step_removelig.execute()
        out_path = os.path.join(self._test_dir, "1uyd_auto.pdb")
        step_removelig.write_generic_by_extension(
            self._test_dir,
            _SGE.PROTEIN_PDB,
        )

        out_file = step_removelig.data.generic.get_files_by_extension("pdb")[
            0
        ].get_data()
        with open(out_path, "w") as f:
            f.write(out_file)
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 126300)
