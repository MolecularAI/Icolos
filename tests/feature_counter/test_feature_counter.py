import unittest

from icolos.core.containers.compound import Compound, Enumeration
from icolos.core.workflow_steps.calculation.feature_counter import StepFeatureCounter
from icolos.utils.enums.program_parameters import FeatureCounterEnum

from icolos.utils.enums.step_enums import StepBaseEnum, StepFeatureCounterEnum

from tests.tests_paths import PATHS_EXAMPLEDATA, get_mol_as_Conformer

_SBE = StepBaseEnum
_FC = FeatureCounterEnum()
_SFC = StepFeatureCounterEnum()


class Test_FeatureCounter(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass

    def setUp(self):
        comp0 = Compound(compound_number=0)
        comp1 = Compound(compound_number=1)
        comp0.add_enumeration(Enumeration(), auto_update=True)
        comp1.add_enumeration(Enumeration(), auto_update=True)
        comp0[0].add_conformers(
            get_mol_as_Conformer(PATHS_EXAMPLEDATA.CLUSTERING_11CONFS), auto_update=True
        )
        comp1[0].add_conformers(
            get_mol_as_Conformer(PATHS_EXAMPLEDATA.SMALL_MOLECULES_SDF_PATH),
            auto_update=True,
        )
        self.comp0 = comp0
        self.comp1 = comp1

    @classmethod
    def tearDownClass(cls):
        pass

    def test_ring_counting(self):
        step_conf = {
            _SBE.STEPID: "01_feature_counting",
            _SBE.STEP_TYPE: _SBE.STEP_FEATURE_COUNTER,
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                },
                _SBE.SETTINGS_ADDITIONAL: {_SFC.FEATURE: _FC.PROPERTY_NUM_RINGS},
            },
        }

        fc_step = StepFeatureCounter(**step_conf)
        fc_step.data.compounds = [self.comp0, self.comp1]

        fc_step.execute()

        self.assertEqual(
            fc_step.get_compounds()[0][0][0]
            .get_molecule()
            .GetProp(_FC.PROPERTY_NUM_RINGS),
            "1",
        )
        self.assertEqual(
            fc_step.get_compounds()[1][0][1]
            .get_molecule()
            .GetProp(_FC.PROPERTY_NUM_RINGS),
            "1",
        )

    def test_aromatic_ring_counting(self):
        step_conf = {
            _SBE.STEPID: "01_feature_counting",
            _SBE.STEP_TYPE: _SBE.STEP_FEATURE_COUNTER,
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                },
                _SBE.SETTINGS_ADDITIONAL: {
                    _SFC.FEATURE: _FC.PROPERTY_NUM_AROMATIC_RINGS
                },
            },
        }

        fc_step = StepFeatureCounter(**step_conf)
        fc_step.data.compounds = [self.comp0, self.comp1]

        fc_step.execute()

        self.assertEqual(
            fc_step.get_compounds()[0][0][0]
            .get_molecule()
            .GetProp(_FC.PROPERTY_NUM_AROMATIC_RINGS),
            "1",
        )
        self.assertEqual(
            fc_step.get_compounds()[1][0][1]
            .get_molecule()
            .GetProp(_FC.PROPERTY_NUM_AROMATIC_RINGS),
            "1",
        )
