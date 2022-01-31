from icolos.core.containers.generic import GenericData
from icolos.utils.enums.program_parameters import ShaepEnum
from tests.tests_paths import (
    PATHS_1UYD,
    PATHS_EXAMPLEDATA,
    get_mol_as_Compound,
    get_mol_as_Conformer,
    MAIN_CONFIG,
)
import unittest
import os

from icolos.utils.enums.step_enums import StepBaseEnum, StepShaepEnum
from icolos.core.workflow_steps.calculation.shaep import StepShaep
from icolos.utils.general.files_paths import attach_root_path

_SBE = StepBaseEnum
_SSE = StepShaepEnum()
_SE = ShaepEnum()


class Test_Shaep(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/shaep")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    def setUp(self):
        # TODO: update to load at least 3 compounds docked (at least 5 poses each)
        mol = get_mol_as_Compound(PATHS_1UYD.NATIVE_LIGAND_SDF)
        conf = get_mol_as_Conformer(PATHS_1UYD.NATIVE_LIGAND_SDF)
        mol[0].add_conformers(conf, auto_update=True)
        self.mol = mol

        with open(PATHS_EXAMPLEDATA.PANTHER_NEGATIVE_IMAGE, "r") as f:
            self.negative_image = f.read()

    def test_shaep(self):
        step_conf = {
            _SBE.STEPID: "01_shaep",
            _SBE.STEP_TYPE: _SBE.STEP_SHAEP,
            _SBE.EXEC: {_SBE.EXEC_BINARYLOCATION: MAIN_CONFIG["SHAEP_LOCATION"]},
        }
        shaep_step = StepShaep(**step_conf)
        shaep_step.data.compounds = [self.mol]
        shaep_step.data.generic.add_file(
            GenericData(file_name="neg_image.mol2", file_data=self.negative_image)
        )
        shaep_step.execute()

        self.assertEqual(
            float(
                shaep_step.get_compounds()[0][0][0]
                .get_molecule()
                .GetProp(_SE.TAG_SHAPE_SIMILARITY)
            ),
            0.663072,
        )
        self.assertEqual(
            float(
                shaep_step.get_compounds()[0][0][0]
                .get_molecule()
                .GetProp(_SE.TAG_ESP_SIMILARITY)
            ),
            0.302026,
        )

        # check, whether the tags got added
        out_path = os.path.join(self._test_dir, "mols_nibr.sdf")
        shaep_step.write_conformers(out_path)
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 2586)
