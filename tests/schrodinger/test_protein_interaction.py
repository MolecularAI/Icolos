import unittest
from icolos.core.workflow_steps.schrodinger.protein_interaction import (
    StepProteinInteraction,
)

from icolos.utils.enums.step_enums import StepBaseEnum
from icolos.core.workflow_steps.schrodinger.prepwizard import StepPrepwizard
from icolos.core.containers.generic import GenericData
from tests.tests_paths import (
    MAIN_CONFIG,
    PATHS_1UYD,
    PATHS_EXAMPLEDATA,
    get_ligands_as_compounds_with_conformers,
    get_mol_as_Compound,
)
import os
from icolos.utils.enums.step_enums import StepBaseEnum, StepGromacsEnum, StepPrepwizEnum
from tests.tests_paths import PATHS_EXAMPLEDATA
from icolos.utils.general.files_paths import attach_root_path
from icolos.core.workflow_steps.schrodinger.prepwizard import StepPrepwizard

_SGE = StepGromacsEnum()
_SBE = StepBaseEnum
_SPE = StepPrepwizEnum()


class TestProteinInteraction(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/protein_interaction")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    def setUp(self):
        base_comp = get_ligands_as_compounds_with_conformers(
            PATHS_EXAMPLEDATA.PROTIEN_INTERACTION_9MER
        )[0]
        self.compounds = [base_comp for _ in range(5)]

    def test_protein_interaction(self):
        step_conf = {
            _SBE.STEPID: "01_protien_interaction",
            _SBE.STEP_TYPE: _SBE.STEP_PROTEIN_INTERACTIONS,
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load schrodinger/2021-4-js-aws",
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {"-group1": "A", "-group2": "A"}
                }
            },
        }

        step_protein_interaction = StepProteinInteraction(**step_conf)
        step_protein_interaction.data.compounds = self.compounds
        step_protein_interaction.execute()
        # extract the dataframe
        out_file = (
            step_protein_interaction.data.compounds[0]
            .get_enumerations()[0]
            .get_conformers()[0]
            .get_extra_data()["interaction_summary"]
        )
        self.assertEqual()
        # out_path = os.path.join(self._test_dir, "test_out.pdb")
        # with open(out_path, "w") as f:
        #     f.write(out_file)
        # stat_inf = os.stat(out_path)
        # self.assertGreater(stat_inf.st_size, 155800)
