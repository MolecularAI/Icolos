from icolos.core.containers.generic import GenericData
import unittest
from icolos.core.workflow_steps.structure_prediction.peptide_embedder import (
    StepPeptideEmbedder,
)
from icolos.utils.general.files_paths import attach_root_path
import os
from tests.tests_paths import PATHS_EXAMPLEDATA
from icolos.utils.enums.step_enums import StepBaseEnum

_SBE = StepBaseEnum


class TestPeptideEmbedder(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/structure_prediction")

        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    def setUp(self):
        with open(attach_root_path(PATHS_EXAMPLEDATA.TEST_FASTA_FILE), "r") as f:
            self.fasta = f.read()

    def test_peptide_embedder(self):
        step_conf = {
            _SBE.STEPID: "01_peptide_embedder",
            _SBE.STEP_TYPE: _SBE.STEP_PEPTIDE_EMBEDDER,
            _SBE.EXEC: {},
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                        "-prime_opt": "OPLS_VERSION=OPLS3e",
                        "-HOST": "cpu-only",
                    },
                },
                _SBE.SETTINGS_ADDITIONAL: {},
            },
        }

        step_embedder = StepPeptideEmbedder(**step_conf)
        fasta_obj = GenericData(file_name="test_seq.fasta", file_data=self.fasta)
        step_embedder.data.generic.add_file(fasta_obj)
        step_embedder.execute()

        out_path = os.path.join(self._test_dir, "sequence_0.pdb")
        step_embedder.write_generic_by_extension(path=self._test_dir, ext="pdb")
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 17504)
