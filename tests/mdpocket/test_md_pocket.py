from icolos.core.containers.generic import GenericData
import unittest
from icolos.utils.enums.step_enums import StepCavExploreEnum, StepBaseEnum
from tests.tests_paths import PATHS_EXAMPLEDATA
from icolos.core.workflow_steps.fpocket.mdpocket import StepMDpocket
from icolos.utils.general.files_paths import attach_root_path

import os

_SBE = StepBaseEnum
_SFP = StepCavExploreEnum()


class Test_MDPocket(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/cavity_explorer")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    def setUp(self):
        # read in the example files here
        # self.cav_folder = PATHS_EXAMPLEDATA.CAVITY_TRJ_FOLDER
        # with open(PATHS_EXAMPLEDATA.CAVITY_DTR_FILE, "rb") as f:
        #     self.dtr_data = f.read()
        # with open(PATHS_EXAMPLEDATA.MD_POCKET_DESMOND_TOP, "r") as f:
        #     self.desmond_pdb = f.read()
        with open(PATHS_EXAMPLEDATA.MDPOCKET_XTC_FILE, "rb") as f:
            self.xtc_data = f.read()
        with open(PATHS_EXAMPLEDATA.MDPOCKET_PDB_FILE, "r") as f:
            self.pdb_file = f.read()

    @classmethod
    def tearDownClass(cls):
        pass

    # def test_MDpocket_desmond(self):
    #     step_conf = {
    #         _SBE.STEPID: "01_cavity_explorer_file_preparation",
    #         _SBE.STEP_TYPE: _SBE.STEP_MDPOCKET,
    #         _SBE.EXEC: {
    #             _SBE.EXEC_PREFIXEXECUTION: "module load fpocket",
    #             _SBE.EXEC_PARALLELIZATION: {
    #                 _SBE.EXEC_PARALLELIZATION_CORES: 8,
    #                 _SBE.EXEC_PARALLELIZATION_MAXLENSUBLIST: 1,
    #             },
    #         },
    #         _SBE.SETTINGS: {
    #             _SBE.SETTINGS_ADDITIONAL: {
    #                 _SFP.SELECTION_TEXT: _SFP.PROTEIN,
    #                 _SFP.TRAJ_TYPE: "desmond",
    #             }
    #         },
    #     }

    #     mdpocket_step = StepMDpocket(**step_conf)
    #     mdpocket_step.data.generic.add_file(
    #         GenericData(
    #             file_name="trj_folder", file_data=self.cav_folder, argument=False
    #         )
    #     )
    #     mdpocket_step.data.generic.add_file(
    #         GenericData(
    #             file_name="structure.pdb", file_data=self.desmond_pdb, argument=True
    #         )
    #     )
    #     mdpocket_step.data.generic.add_file(
    #         GenericData(file_name="clickme.dtr", file_data=self.dtr_data, argument=True)
    #     )
    #     mdpocket_step.execute()

    #     out_path = os.path.join(self._test_dir, "pocket_0_descriptors.txt")
    #     mdpocket_step.write_generic_by_extension(self._test_dir, "txt")

    #     stat_inf = os.stat(out_path)
    #     self.assertGreater(stat_inf.st_size, 187400)

    def test_MDpocket_xtc(self):
        step_conf = {
            _SBE.STEPID: "01_cavity_explorer_file_preparation",
            _SBE.STEP_TYPE: _SBE.STEP_MDPOCKET,
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load fpocket",
                _SBE.EXEC_PARALLELIZATION: {
                    _SBE.EXEC_PARALLELIZATION_CORES: 4,
                    _SBE.EXEC_PARALLELIZATION_MAXLENSUBLIST: 2,  # this will be automatically overwritten
                },
                _SBE.EXEC_FAILUREPOLICY: {_SBE.EXEC_FAILUREPOLICY_NTRIES: 1},
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ADDITIONAL: {
                    _SFP.SELECTION_TEXT: _SFP.PROTEIN,
                    _SFP.TRAJ_TYPE: "gromacs",
                }
            },
        }

        mdpocket_step = StepMDpocket(**step_conf)

        mdpocket_step.data.generic.add_file(
            GenericData(
                file_name="structure.xtc", file_data=self.xtc_data, argument=True
            )
        )
        mdpocket_step.data.generic.add_file(
            GenericData(
                file_name="structure.pdb", file_data=self.pdb_file, argument=True
            )
        )
        mdpocket_step.execute()

        out_path = os.path.join(self._test_dir, "pocket_1_descriptors.txt")
        mdpocket_step.write_generic_by_name(self._test_dir, "pocket_1_descriptors.txt")
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 700)
