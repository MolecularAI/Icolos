import unittest

from icolos.core.containers.compound import Compound, Enumeration
from icolos.core.workflow_steps.calculation.rms_filter import StepRMSFilter
from icolos.utils.enums.step_enums import StepBaseEnum, StepRMSFilterEnum
from tests.tests_paths import PATHS_EXAMPLEDATA, get_mol_as_Conformer

_SBE = StepBaseEnum
_SRF = StepRMSFilterEnum()


class Test_RMSfilter(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass

    def setUp(self):
        pass

    @classmethod
    def tearDownClass(cls):
        pass

    def test_RMSfiltering_alignmol_descending(self):
        step_conf = {
            _SBE.STEPID: "01_RMSfiltering",
            _SBE.STEP_TYPE: _SBE.STEP_RMSFILTER,
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                },
                _SBE.SETTINGS_ADDITIONAL: {
                    _SRF.METHOD: _SRF.METHOD_ALIGNMOL,
                    _SRF.THRESHOLD: 1,
                    _SRF.ORDER_BY: "E_cosmo",
                    _SRF.ORDER_ASCENDING: False,
                },
            },
        }

        rf_step = StepRMSFilter(**step_conf)
        rf_step.get_compounds().append(Compound(compound_number=1))
        rf_step.get_compounds()[0].add_enumeration(Enumeration(), auto_update=True)
        conformers = get_mol_as_Conformer(PATHS_EXAMPLEDATA.CLUSTERING_11CONFS)
        rf_step.data.compounds[0][0].add_conformers(conformers, auto_update=True)

        self.assertEqual(len(rf_step.get_compounds()[0][0].get_conformers()), 11)
        rf_step.execute()
        self.assertEqual(len(rf_step.get_compounds()[0][0].get_conformers()), 7)

        self.assertListEqual(
            list(
                float(
                    rf_step.get_compounds()[0][0]
                    .get_conformers()[i]
                    .get_molecule()
                    .GetProp("E_cosmo")
                )
                for i in range(7)
            ),
            [
                -943306.7731,
                -943304.5548,
                -943301.0009,
                -943300.9934,
                -943303.7802,
                -943304.0485,
                -943304.0517,
            ],
        )

        step_conf[_SBE.SETTINGS][_SBE.SETTINGS_ADDITIONAL][_SRF.THRESHOLD] = 1.5
        rf_step = StepRMSFilter(**step_conf)
        rf_step.get_compounds().append(Compound(compound_number=1))
        rf_step.get_compounds()[0].add_enumeration(Enumeration(), auto_update=True)
        conformers = get_mol_as_Conformer(PATHS_EXAMPLEDATA.CLUSTERING_11CONFS)
        rf_step.data.compounds[0][0].add_conformers(conformers, auto_update=True)

        self.assertEqual(len(rf_step.get_compounds()[0][0].get_conformers()), 11)
        rf_step.execute()
        self.assertEqual(len(rf_step.get_compounds()[0][0].get_conformers()), 3)

    def test_RMSfiltering_alignmol_ascending(self):
        step_conf = {
            _SBE.STEPID: "01_RMSfiltering",
            _SBE.STEP_TYPE: _SBE.STEP_RMSFILTER,
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                },
                _SBE.SETTINGS_ADDITIONAL: {
                    _SRF.METHOD: _SRF.METHOD_ALIGNMOL,
                    _SRF.THRESHOLD: 1,
                    _SRF.ORDER_BY: "E_cosmo",
                    _SRF.ORDER_ASCENDING: True,
                },
            },
        }

        rf_step = StepRMSFilter(**step_conf)
        rf_step.get_compounds().append(Compound(compound_number=1))
        rf_step.get_compounds()[0].add_enumeration(Enumeration(), auto_update=True)
        conformers = get_mol_as_Conformer(PATHS_EXAMPLEDATA.CLUSTERING_11CONFS)
        rf_step.data.compounds[0][0].add_conformers(conformers, auto_update=True)

        self.assertEqual(len(rf_step.get_compounds()[0][0].get_conformers()), 11)
        rf_step.execute()
        self.assertEqual(len(rf_step.get_compounds()[0][0].get_conformers()), 5)

        self.assertListEqual(
            list(
                float(
                    rf_step.get_compounds()[0][0]
                    .get_conformers()[i]
                    .get_molecule()
                    .GetProp("E_cosmo")
                )
                for i in range(5)
            ),
            [-943304.5487, -943300.2823, -943303.7733, -943304.0485, -943304.0517],
        )

    def test_RMSfiltering_best(self):
        step_conf = {
            _SBE.STEPID: "01_RMSfiltering",
            _SBE.STEP_TYPE: _SBE.STEP_RMSFILTER,
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                },
                _SBE.SETTINGS_ADDITIONAL: {
                    _SRF.METHOD: _SRF.METHOD_BEST,
                    _SRF.THRESHOLD: 1,
                    _SRF.ORDER_BY: "E_cosmo",
                    _SRF.ORDER_ASCENDING: False,
                },
            },
        }

        rf_step = StepRMSFilter(**step_conf)
        rf_step.get_compounds().append(Compound(compound_number=1))
        rf_step.get_compounds()[0].add_enumeration(Enumeration(), auto_update=True)
        conformers = get_mol_as_Conformer(PATHS_EXAMPLEDATA.CLUSTERING_11CONFS)
        rf_step.data.compounds[0][0].add_conformers(conformers, auto_update=True)

        self.assertEqual(len(rf_step.get_compounds()[0][0].get_conformers()), 11)
        rf_step.execute()
        self.assertEqual(len(rf_step.get_compounds()[0][0].get_conformers()), 4)

        self.assertListEqual(
            list(
                float(
                    rf_step.get_compounds()[0][0]
                    .get_conformers()[i]
                    .get_molecule()
                    .GetProp("E_cosmo")
                )
                for i in range(4)
            ),
            [-943306.7731, -943304.5548, -943301.0009, -943304.0517],
        )

    def test_RMSfiltering_best_notordered(self):
        step_conf = {
            _SBE.STEPID: "01_RMSfiltering",
            _SBE.STEP_TYPE: _SBE.STEP_RMSFILTER,
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                },
                _SBE.SETTINGS_ADDITIONAL: {
                    _SRF.METHOD: _SRF.METHOD_BEST,
                    _SRF.THRESHOLD: 1,
                },
            },
        }

        rf_step = StepRMSFilter(**step_conf)
        rf_step.get_compounds().append(Compound(compound_number=1))
        rf_step.get_compounds()[0].add_enumeration(Enumeration(), auto_update=True)
        conformers = get_mol_as_Conformer(PATHS_EXAMPLEDATA.CLUSTERING_11CONFS)
        rf_step.data.compounds[0][0].add_conformers(conformers, auto_update=True)

        self.assertEqual(len(rf_step.get_compounds()[0][0].get_conformers()), 11)
        rf_step.execute()
        self.assertEqual(len(rf_step.get_compounds()[0][0].get_conformers()), 4)

        self.assertListEqual(
            list(
                float(
                    rf_step.get_compounds()[0][0]
                    .get_conformers()[i]
                    .get_molecule()
                    .GetProp("E_cosmo")
                )
                for i in range(4)
            ),
            [-943306.7731, -943304.5487, -943301.0009, -943304.0485],
        )
