import unittest
import os

from icolos.core.containers.compound import Compound, Enumeration
from icolos.core.workflow_steps.calculation.boltzmann_weighting import (
    StepBoltzmannWeighting,
)

from icolos.utils.enums.step_enums import StepBaseEnum, StepBoltzmannWeightingEnum

from tests.tests_paths import PATHS_EXAMPLEDATA, get_mol_as_Conformer
from icolos.utils.general.files_paths import attach_root_path

_SBE = StepBaseEnum
_SBWE = StepBoltzmannWeightingEnum()


class Test_BoltzmannWeighting(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/BoltzmannWeighting")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    def setUp(self):
        # this file has the necessary properties for the different solvents annotated as tags
        self._example_mol_path = (
            PATHS_EXAMPLEDATA.EPSA_BOLTZMANN_WEIGHTING_EXAMPLE_MOLECULE
        )

    @classmethod
    def tearDownClass(cls):
        pass

    def test_predict_ePSA_with_descriptors(self):
        step_conf = {
            _SBE.STEPID: "01_boltzmann_weighting",
            _SBE.STEP_TYPE: _SBE.STEP_BOLTZMANN_WEIGHTING,
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                        _SBWE.PROPERTIES: [
                            {
                                _SBWE.PROPERTIES_INPUT: "G_h2o",
                                _SBWE.PROPERTIES_OUTPUT: "boltzfactor_wat",
                            },
                            {
                                _SBWE.PROPERTIES_INPUT: "G_meoh",
                                _SBWE.PROPERTIES_OUTPUT: "boltzfactor_meoh",
                            },
                            {
                                _SBWE.PROPERTIES_INPUT: "G_octanol",
                                _SBWE.PROPERTIES_OUTPUT: "boltzfactor_octanol",
                            },
                            {
                                _SBWE.PROPERTIES_INPUT: "G_dmso",
                                _SBWE.PROPERTIES_OUTPUT: "boltzfactor_dmso",
                            },
                            {
                                _SBWE.PROPERTIES_INPUT: "G_cychex",
                                _SBWE.PROPERTIES_OUTPUT: "boltzfactor_cychex",
                            },
                            {
                                _SBWE.PROPERTIES_INPUT: "G_chcl3",
                                _SBWE.PROPERTIES_OUTPUT: "boltzfactor_chcl3",
                            },
                            {
                                _SBWE.PROPERTIES_INPUT: "G_acn",
                                _SBWE.PROPERTIES_OUTPUT: "boltzfactor_acn",
                            },
                            {
                                _SBWE.PROPERTIES_INPUT: "G_thf",
                                _SBWE.PROPERTIES_OUTPUT: "boltzfactor_thf",
                            },
                        ],
                        _SBWE.WEIGHT: {
                            _SBWE.WEIGHT_INPUT: [
                                "area",
                                "HB_acc",
                                "volume",
                                "HB_don",
                                "sigma2",
                                "Gsolv_meoh",
                            ],
                            _SBWE.WEIGHT_OUTPUT_PREFIX: "bf_weighted",
                            _SBWE.WEIGHT_PROPERTIES: [
                                "boltzfactor_dmso",
                                "boltzfactor_wat",
                                "boltzfactor_meoh",
                                "boltzfactor_cychex",
                            ],
                        },
                    }
                }
            },
        }
        bweigh_step = StepBoltzmannWeighting(**step_conf)
        bweigh_step.get_compounds().append(Compound())
        bweigh_step.get_compounds()[0].add_enumeration(Enumeration(), auto_update=True)
        conformer = get_mol_as_Conformer(self._example_mol_path)
        bweigh_step.data.compounds[0][0].add_conformers(conformer, auto_update=True)
        bweigh_step.execute()

        self.assertEqual(len(bweigh_step.get_compounds()), 1)
        self.assertEqual(len(bweigh_step.get_compounds()[0]), 1)
        self.assertEqual(len(bweigh_step.get_compounds()[0][0]), 1)

        self.assertListEqual(
            list(
                bweigh_step.get_compounds()[0][0]
                .get_conformers()[0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [1.9524, -0.9976, -1.5113],
        )
        self.assertEqual(
            19,
            bweigh_step.get_compounds()[0][0]
            .get_conformers()[0]
            .get_molecule()
            .GetNumAtoms(),
        )

        # check SDF write-out (including Boltzmann-weighted properties as tags)
        out_path = os.path.join(self._test_dir, "boltzmann_weighted_annotated.sdf")
        bweigh_step.write_conformers(out_path)
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 4419)
