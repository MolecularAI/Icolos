from icolos.core.containers.generic import GenericData
import unittest
import os
from icolos.core.workflow_steps.schrodinger.fep_analysis import StepFepPlusAnalysis
from icolos.utils.enums.step_enums import StepBaseEnum, StepFepPlusEnum, StepGlideEnum
from tests.tests_paths import (
    PATHS_EXAMPLEDATA,
    get_ligands_as_compounds_with_conformers,
    PATHS_1UYD,
)
from icolos.utils.general.files_paths import attach_root_path

_SBE = StepBaseEnum
_SFE = StepFepPlusEnum()
_SGE = StepGlideEnum()


class Test_FepPlusAnalysis(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/fep_plus")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    def setUp(self):
        with open(attach_root_path(PATHS_EXAMPLEDATA.FEP_PLUS_MULTISIM_LONG), "r") as f:
            self.log = f.read()

        self.mol = get_ligands_as_compounds_with_conformers(
            attach_root_path(PATHS_1UYD.LIG_SDF)
        )

    def test_fep_analysis(self):
        step_conf = {
            _SBE.STEPID: "test_fep_analysis",
            _SBE.STEP_TYPE: "fep_analysis",
            _SBE.EXEC: {_SBE.EXEC_PREFIXEXECUTION: "module load schrodinger/2020-4"},
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {},
                _SBE.SETTINGS_ADDITIONAL: {_SFE.REFERENCE_DG: -10.76},
            },
        }

        step_fep_analysis = StepFepPlusAnalysis(**step_conf)
        step_fep_analysis.data.compounds = self.mol
        step_fep_analysis.data.generic.add_file(
            GenericData(
                file_name="test_multisim.log", file_data=self.log, argument=True
            )
        )
        step_fep_analysis.execute()
        # now confirm that the values have been parsed out of the log file properly
        # true_conf_energies = ['2.67+-0.48', '0.00+-0.40', '2.86+-0.42', '8.88+-0.52', '3.09+-0.41']
        true_conf_energies = [
            -10.76,
            -6.72,
            -8.87,
            -7.1,
            -7.36,
            -9.18,
            -10.38,
            -9.2,
            -5.73,
            -7.91,
            -9.16,
            -7.38,
            -7.44,
            -1.92,
            -6.78,
            -6.35,
            -2.54,
            -7.17,
            -6.89,
            -8.32,
            -8.21,
            -6.92,
            -6.28,
            -7.03,
            -8.23,
            -11.38,
            -9.14,
            -7.35,
            -7.21,
            -7.39,
            -1.48,
            -8.02,
            -7.14,
            -6.3,
            -7.59,
            -9.79,
            -6.84,
            -7.1,
        ]
        conformer_energies = []
        for compound in step_fep_analysis.data.compounds:
            conformer_energies.append(
                compound.get_enumerations()[0]
                .get_conformers()[0]
                .get_molecule()
                .GetProp("map_dG")
            )
        for idx, energy in enumerate(conformer_energies):
            self.assertAlmostEqual(
                float(energy.split("+-")[0]), true_conf_energies[idx], 2
            )
