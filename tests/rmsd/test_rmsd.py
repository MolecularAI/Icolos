import unittest
from copy import deepcopy
from typing import List

from rdkit.Geometry.rdGeometry import Point3D

from icolos.core.containers.compound import Compound, Enumeration, unroll_conformers
from icolos.core.workflow_steps.calculation.rmsd import StepRMSD

from icolos.utils.enums.step_enums import (
    StepBaseEnum,
    StepRMSDEnum,
    StepDataManipulationEnum,
)

from tests.tests_paths import PATHS_EXAMPLEDATA, get_mol_as_Conformer

_SBE = StepBaseEnum
_SR = StepRMSDEnum()
_SDM = StepDataManipulationEnum()


def _match_as_generic(
    comp_list_1: List[Compound], comp_list_2: List[Compound]
) -> List[Compound]:
    comp2_conf_unrolled = unroll_conformers(comp_list_2)

    # attach the second version of the conformers as generic field to the "real" input
    # (as would be done by the data manipulator)
    for comp in comp_list_1:
        for enum in comp:
            for conf in enum:
                conf.add_extra_data(
                    key=_SDM.KEY_MATCHED,
                    data=[
                        c
                        for c in comp2_conf_unrolled
                        if conf.get_index_string() == c.get_index_string()
                    ],
                )
    return comp_list_1


class Test_RMSD(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass

    def setUp(self):
        conformers = get_mol_as_Conformer(PATHS_EXAMPLEDATA.CLUSTERING_11CONFS)

        # Compound 1 with 1 enumeration and 11 conformers
        self.comp1 = Compound(compound_number=1)
        self.comp1.add_enumeration(Enumeration(), auto_update=True)
        self.comp1[0].add_conformers(deepcopy(conformers), auto_update=True)

        # Compound 2 with 1 enumeration and 11 conformers, change of some coordinates
        self.comp2 = Compound(compound_number=1)
        self.comp2.add_enumeration(Enumeration(), auto_update=True)
        self.comp2[0].add_conformers(deepcopy(conformers), auto_update=True)
        self.comp2[0][1].get_molecule().GetConformer().SetAtomPosition(
            0, Point3D(-4.2239, -0.441, 0.2458)
        )
        self.comp2[0][7].get_molecule().GetConformer().SetAtomPosition(
            0, Point3D(-1.5442, -0.7854, 0.5883)
        )

    @classmethod
    def tearDownClass(cls):
        pass

    def test_RMSD_conformers_matched(self):
        step_conf = {
            _SBE.STEPID: "01_RMSD",
            _SBE.STEP_TYPE: _SBE.STEP_RMSD,
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                },
                _SBE.SETTINGS_ADDITIONAL: {_SR.METHOD: _SR.METHOD_ALIGNMOL},
            },
        }

        rf_step = StepRMSD(**step_conf)
        rf_step.get_compounds().append(_match_as_generic([self.comp1], [self.comp2])[0])
        self.assertEqual(len(rf_step.get_compounds()[0][0][0].get_extra_data()), 1)

        rf_step.execute()

        self.assertEqual(
            rf_step.get_compounds()[0][0][1].get_molecule().GetProp(_SR.RMSD_TAG),
            "8.268",
        )
        self.assertEqual(
            rf_step.get_compounds()[0][0][1]
            .get_extra_data()[_SDM.KEY_MATCHED][0]
            .get_molecule()
            .GetProp(_SR.RMSD_TAG),
            "8.268",
        )
        self.assertEqual(
            rf_step.get_compounds()[0][0][3].get_molecule().GetProp(_SR.RMSD_TAG), "0.0"
        )
        self.assertEqual(
            rf_step.get_compounds()[0][0][3]
            .get_extra_data()[_SDM.KEY_MATCHED][0]
            .get_molecule()
            .GetProp(_SR.RMSD_TAG),
            "0.0",
        )
        self.assertEqual(len(rf_step.get_compounds()[0][0].get_conformers()), 11)
