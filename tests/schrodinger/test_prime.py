import unittest
import os

from icolos.core.workflow_steps.schrodinger.prime import StepPrime

from icolos.utils.enums.step_enums import StepBaseEnum, StepPrimeEnum, TokenGuardEnum
from icolos.utils.enums.program_parameters import PrimeEnum

from tests.tests_paths import (
    PATHS_EXAMPLEDATA,
    get_mol_as_Compound,
    get_ligands_as_compounds_with_conformers,
)
from icolos.utils.general.files_paths import attach_root_path

_SBE = StepBaseEnum
_SPE = StepPrimeEnum()
_CE = PrimeEnum()
_TE = TokenGuardEnum()


class Test_Prime(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/prime_test")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    def setUp(self):
        with open(attach_root_path(PATHS_EXAMPLEDATA.PRIME_POSEVIEWER), "rb") as f:
            self._poseviewer = f.read()
        self._molecule = get_mol_as_Compound(
            attach_root_path(PATHS_EXAMPLEDATA.PRIME_DOCKED_LIGAND_SDF)
        )
        self._conformers = get_ligands_as_compounds_with_conformers(
            attach_root_path(PATHS_EXAMPLEDATA.LIGANDS_1UYD)
        )

    @classmethod
    def tearDownClass(cls):
        pass

    def test_Prime_run(self):
        # TODO: make sure the original execution mode (on enumerations) works ok
        # * Pull the molecule from the enumeration if no conformers attached
        # * add conformer to the enum at the end
        step_conf = {
            _SBE.STEPID: "01_prime",
            _SBE.STEP_TYPE: _SBE.STEP_PRIME,
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load schrodinger/2021-2-js-aws",
                _SBE.EXEC_PARALLELIZATION: {
                    _SBE.EXEC_PARALLELIZATION_CORES: 4,
                    _SBE.EXEC_PARALLELIZATION_MAXLENSUBLIST: 2,
                },
                _SBE.EXEC_FAILUREPOLICY: {_SBE.EXEC_FAILUREPOLICY_NTRIES: 1},
            },
            _TE.TG: {
                _TE.TG_PREFIX_EXECUTION: "module load schrodinger/2021-2-js-aws",
                _TE.TG_TOKEN_POOLS: {"PRIMEX_MAIN": 8},
                _TE.TG_WAIT_INTERVAL_SECONDS: 30,
                _TE.TG_WAIT_LIMIT_SECONDS: 900,
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                        "-prime_opt": "OPLS_VERSION=OPLS3e",
                        "-HOST": "cpu-only",
                    },
                },
                _SBE.SETTINGS_ADDITIONAL: {
                    _SPE.RECEPTOR: attach_root_path(PATHS_EXAMPLEDATA.RECEPTOR_1UYD)
                },
            },
        }

        prime_step = StepPrime(**step_conf)
        prime_step.data.compounds = [self._molecule]
        prime_step.execute()

        self.assertEqual(len(prime_step.get_compounds()[0][0].get_conformers()), 1)
        # molecule coordinates should not be touched by the execution (conformer is optimized though)
        self.assertListEqual(
            list(
                prime_step.get_compounds()[0][0]
                .get_molecule()
                .GetConformer(0)
                .GetPositions()[0]
            ),
            [15.2886, 52.7, 69.7128],
        )

        # check write-out
        out_path = os.path.join(self._test_dir, "prime_output_file.sdf")
        prime_step.write_conformers(out_path)
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 10000)
        self.assertGreater(13500, stat_inf.st_size)

    def test_prime_run_conformers(self):
        step_conf = {
            _SBE.STEPID: "01_prime",
            _SBE.STEP_TYPE: _SBE.STEP_PRIME,
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load schrodinger/2021-2-js-aws",
                _SBE.EXEC_PARALLELIZATION: {
                    _SBE.EXEC_PARALLELIZATION_CORES: 32,
                    _SBE.EXEC_PARALLELIZATION_MAXLENSUBLIST: 1,
                },
                _SBE.EXEC_FAILUREPOLICY: {_SBE.EXEC_FAILUREPOLICY_NTRIES: 1},
            },
            _TE.TG: {
                _TE.TG_PREFIX_EXECUTION: "module load schrodinger/2021-2-js-aws",
                _TE.TG_TOKEN_POOLS: {"PRIMEX_MAIN": 8},
                _TE.TG_WAIT_INTERVAL_SECONDS: 30,
                _TE.TG_WAIT_LIMIT_SECONDS: 900,
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                        "-prime_opt": "OPLS_VERSION=OPLS3e",
                        "-HOST": "cpu-only",
                    },
                },
                _SBE.SETTINGS_ADDITIONAL: {
                    _SPE.RECEPTOR: attach_root_path(PATHS_EXAMPLEDATA.RECEPTOR_1UYD)
                },
            },
        }
        prime_step = StepPrime(**step_conf)
        prime_step.data.compounds = self._conformers
        prime_step.execute()
        out_path = os.path.join(self._test_dir, "prime_conformers_output.sdf")
        prime_step.write_conformers(out_path)
        scores = [
            "-46.4912412523436",
            "-49.5744863214668",
            "-63.4520243626994",
            "-55.2546247037599",
            "-35.0457131568983",
            "-37.584671678831",
            "-52.3315306739823",
            "-42.1457765778323",
            "-39.0962071597705",
            "-46.9267618228951",
            "-41.4015029031088",
            "-49.0027294452047",
            "-45.297078493255",
            "-47.1669750502297",
            "-50.2110899116497",
            "-38.8494636817877",
            "-41.6326792228592",
            "-43.6924482130898",
            "-46.738882435201",
            "-45.242419676907",
            "-36.5693940219298",
            "-57.9606138506851",
            "-55.4918326231546",
            "-39.724716804717",
            "-50.0105377772616",
            "-46.9162249942074",
            "-46.2790546176639",
            "-43.8232309398354",
            "-49.7540870967205",
            "-53.7133446915177",
            "-51.6633994627191",
            "-54.2858218610409",
            "-42.9129639283819",
            "-49.1980564160085",
            "-52.7421500005312",
            "-50.953927771995",
            "-59.8079546364734",
            "-53.20869108637",
            "-42.9971732771755",
            "-46.3393621442165",
            "-39.1124509414121",
            "-26.9291589283248",
            "-48.0546634882376",
            "-58.0973312599281",
            "-52.8690868697358",
        ]
        flattened_conformers_scores = []
        for compound in prime_step.data.compounds:
            for enumeration in compound.get_enumerations():
                for conformer in enumeration.get_conformers():
                    flattened_conformers_scores.append(
                        conformer.get_molecule().GetProp(_SPE.MMGBSA_SCORE)
                    )
        # self.assertEqual(float(prime_step.get_compounds()[0].get_enumerations()[0].get_conformers()[0].get_molecule()\
        #                  .GetProp('r_psp_MMGBSA_dG_Bind')), -69.9651350867098)

        for trial, value in zip(flattened_conformers_scores, scores):
            self.assertEqual(round(float(trial)), round(float(value)))
