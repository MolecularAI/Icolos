import os
import unittest
from icolos.core.composite_agents.workflow import WorkFlow
from icolos.core.step_utils.input_preparator import (
    InputPreparator,
    StepInputParameters,
    StepInputSource,
)
from icolos.core.containers.compound import Conformer, Enumeration, Compound
from icolos.core.workflow_steps.step import StepBase
from icolos.utils.enums.step_enums import StepBaseEnum
from icolos.utils.general.files_paths import attach_root_path
from tests.tests_paths import PATHS_EXAMPLEDATA

_SBE = StepBaseEnum


class Test_InputPreparator(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/InputPreparator")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    def setUp(self):
        # comp1 has 2 enumerations, one with 2 and one with 3 conformers
        comp1 = Compound(name="test_molecule", compound_number=0)
        comp1_enum1 = Enumeration(smile="abc", molecule=None, enumeration_id=1)
        comp1_enum1.add_conformer(Conformer(conformer_id=0), auto_update=True)
        comp1_enum1.add_conformer(Conformer(conformer_id=0), auto_update=True)
        comp1_enum2 = Enumeration(smile="def", molecule=None, enumeration_id=2)
        comp1_enum2.add_conformer(Conformer(conformer_id=0), auto_update=True)
        comp1_enum2.add_conformer(Conformer(conformer_id=0), auto_update=True)
        comp1_enum2.add_conformer(Conformer(conformer_id=0), auto_update=True)
        comp1.add_enumeration(comp1_enum1, auto_update=False)
        comp1.add_enumeration(comp1_enum2, auto_update=False)

        source1 = StepInputSource(
            source="mol1:cccccc1",
            source_type=_SBE.INPUT_SOURCE_TYPE_STRING,
            source_field="new_string",
        )
        source2 = StepInputSource(
            source="prev_step", source_type=_SBE.INPUT_SOURCE_TYPE_STEP
        )
        source3 = StepInputSource(
            source="mock_step",
            source_type=_SBE.INPUT_SOURCE_TYPE_STEP,
            source_field="old_input_field",
            target_field="new_input_field",
        )
        source4 = StepInputSource(
            source="mol2:cccc1", source_type=_SBE.INPUT_SOURCE_TYPE_STRING
        )
        source5 = StepInputSource(
            source=attach_root_path(PATHS_EXAMPLEDATA.PARACETAMOL_COSMO),
            source_type=_SBE.INPUT_SOURCE_TYPE_PATH,
            source_field="cosmo",
            target_field="cosmo",
        )
        source6 = StepInputSource(
            source=attach_root_path(PATHS_EXAMPLEDATA.PARACETAMOL_COSMO),
            source_type=_SBE.INPUT_SOURCE_TYPE_FILE,
            source_field="cosmo_filepath",
            target_field="cosmo_test_file",
        )
        source7 = StepInputSource(
            source=PATHS_EXAMPLEDATA.PANTHER_NEGATIVE_IMAGE,
            extension="mol2",
        )
        self.params = StepInputParameters(
            compounds=[source1, source4, source2], generic=[source7]
        )
        blank_params = StepInputParameters(compounds=[], generic=[])
        mock_step = StepBase(step_id="mock_step", type=None, input=self.params)
        prev_step = StepBase(step_id="prev_step", type=None, input=blank_params)
        prev_step.data.compounds = [comp1]

        workflow = WorkFlow()
        workflow.add_step(prev_step)
        workflow.add_step(mock_step)
        self.workflow = workflow

    @classmethod
    def tearDownClass(cls):
        pass

    def test_input_preparation(self):
        preparator = InputPreparator(workflow=self.workflow, logger=None)
        data, work_dir = preparator.generate_input(
            step_input=self.params, step_type=_SBE.STEP_SHAEP
        )
        self.assertEqual(len(data.compounds), 3)
        self.assertEqual(len(data.generic.get_all_files()), 1)
        print(data.generic.get_all_files())
        with open(PATHS_EXAMPLEDATA.PANTHER_NEGATIVE_IMAGE, "r") as f:
            file = f.read()
        self.assertEqual(
            data.generic.get_file_by_name("1uyd_negative_image.mol2").get_data(), file
        )
        self.assertEqual(len(data.compounds[1]), 1)
        self.assertEqual((len(data.compounds[2][1])), 3)
