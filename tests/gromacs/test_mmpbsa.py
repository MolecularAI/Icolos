from icolos.core.workflow_steps.gromacs.mmpbsa import StepGMXmmpbsa
from icolos.core.containers.generic import GenericData
import unittest
import os
from icolos.utils.enums.step_enums import StepBaseEnum, StepGromacsEnum
from tests.tests_paths import PATHS_EXAMPLEDATA, export_unit_test_env_vars, MAIN_CONFIG
from icolos.utils.general.files_paths import attach_root_path
from time import time

_SGE = StepGromacsEnum()
_SBE = StepBaseEnum


class Test_MMPBSA(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls._test_dir = attach_root_path("tests/junk/gromacs")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

        export_unit_test_env_vars()

    def setUp(self) -> None:
        with open(PATHS_EXAMPLEDATA.GROMACS_HOLO_STRUCTURE_GRO, "r") as f:
            self.structure = f.read()

        with open(PATHS_EXAMPLEDATA.MMPBSA_TOP, "r") as f:
            self.topol = f.read()

        with open(PATHS_EXAMPLEDATA.MMPBSA_TPR, "rb") as f:
            self.tpr_file = f.read()

        with open(PATHS_EXAMPLEDATA.MMPBSA_XTC, "rb") as f:
            self.xtc_file = f.read()

        with open(PATHS_EXAMPLEDATA.MMPBSA_POSRE, "rb") as f:
            self.posre = f.read()

        with open(PATHS_EXAMPLEDATA.MMPBSA_LIG_ITP, "rb") as f:
            self.lig_itp = f.read()

        with open(PATHS_EXAMPLEDATA.MMPBSA_LIG_POSRE, "rb") as f:
            self.lig_posre = f.read()

    def test_protein_lig_single_traj(self):
        step_conf = {
            _SBE.STEPID: "test_gmmpbsa",
            _SBE.STEP_TYPE: "gmx_mmpbsa",
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a && module load gmx_MMPBSA && module load AmberTools/21-fosscuda-2019a-Python-3.7.2"
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                },
                _SBE.SETTINGS_ADDITIONAL: {
                    _SGE.FORCEFIELD: MAIN_CONFIG["FORCEFIELD"],
                    _SGE.COUPLING_GROUPS: "Protein Other",
                },
            },
        }
        step_mmpbsa = StepGMXmmpbsa(**step_conf)
        step_mmpbsa.data.generic.add_file(
            GenericData(file_name="structure.gro", file_data=self.structure)
        )
        step_mmpbsa.data.generic.add_file(
            GenericData(file_name="topol.top", file_data=self.topol)
        )
        step_mmpbsa.data.generic.add_file(
            GenericData(file_name="structure.xtc", file_data=self.xtc_file)
        )
        step_mmpbsa.data.generic.add_file(
            GenericData(file_name="structure.tpr", file_data=self.tpr_file)
        )
        step_mmpbsa.data.generic.add_file(
            GenericData(file_name="posre.itp", file_data=self.posre)
        )
        step_mmpbsa.data.generic.add_file(
            GenericData(file_name="DMP:100.itp", file_data=self.lig_itp)
        )
        step_mmpbsa.data.generic.add_file(
            GenericData(file_name="posre_DMP:100.itp", file_data=self.lig_posre)
        )
        step_mmpbsa.execute()
        out_path = os.path.join(self._test_dir, "FINAL_RESULTS_MMPBSA.dat")
        step_mmpbsa.write_generic_by_extension(self._test_dir, "dat")
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 5570)

    def test_protein_lig_single_traj_custom_file(self):

        step_conf = {
            _SBE.STEPID: "test_gmmpbsa",
            _SBE.STEP_TYPE: "gmx_mmpbsa",
            _SBE.EXEC: {
                _SBE.EXEC_PREFIXEXECUTION: "module load GROMACS/2020.3-fosscuda-2019a && module load gmx_MMPBSA && module load AmberTools/21-fosscuda-2019a-Python-3.7.2"
            },
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                },
                _SBE.SETTINGS_ADDITIONAL: {
                    _SGE.FORCEFIELD: MAIN_CONFIG["FORCEFIELD"],
                    _SGE.COUPLING_GROUPS: "Protein Other",
                    _SGE.INPUT_FILE: PATHS_EXAMPLEDATA.MMPBSA_CUSTOM_INPUT,
                    "ntasks": 2,
                },
            },
        }
        step_mmpbsa = StepGMXmmpbsa(**step_conf)
        step_mmpbsa.data.generic.add_file(
            GenericData(file_name="structure.gro", file_data=self.structure)
        )
        step_mmpbsa.data.generic.add_file(
            GenericData(file_name="topol.top", file_data=self.topol)
        )
        step_mmpbsa.data.generic.add_file(
            GenericData(file_name="structure.xtc", file_data=self.xtc_file)
        )
        step_mmpbsa.data.generic.add_file(
            GenericData(file_name="structure.tpr", file_data=self.tpr_file)
        )
        step_mmpbsa.data.generic.add_file(
            GenericData(file_name="posre.itp", file_data=self.posre)
        )
        step_mmpbsa.data.generic.add_file(
            GenericData(file_name="DMP:100.itp", file_data=self.lig_itp)
        )
        step_mmpbsa.data.generic.add_file(
            GenericData(file_name="posre_DMP:100.itp", file_data=self.lig_posre)
        )
        t1 = time()
        step_mmpbsa.execute()
        exec_time = time() - t1
        print("single traj exec time, custom input", exec_time)
        out_path = os.path.join(self._test_dir, "FINAL_RESULTS_MMPBSA.dat")
        step_mmpbsa.write_generic_by_extension(self._test_dir, "dat")
        stat_inf = os.stat(out_path)

        self.assertGreater(stat_inf.st_size, 4680)
