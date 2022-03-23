from icolos.core.containers.gmx_state import GromacsState
from icolos.core.workflow_steps.gromacs.mmpbsa import StepGMXmmpbsa
from icolos.core.containers.generic import GenericData
import unittest
import os
from icolos.utils.enums.step_enums import StepBaseEnum, StepGromacsEnum
from tests.tests_paths import PATHS_EXAMPLEDATA, export_unit_test_env_vars, MAIN_CONFIG
from icolos.utils.general.files_paths import attach_root_path
import numpy as np

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
            self.structure = f.readlines()

        with open(PATHS_EXAMPLEDATA.GROMACS_1BVG_TPR, "rb") as f:
            self.tpr_file = f.read()

        with open(PATHS_EXAMPLEDATA.GROMACS_1BVG_XTC, "rb") as f:
            self.xtc_file = f.read()
        with open(PATHS_EXAMPLEDATA.MMPBSA_POSRE, "rb") as f:
            self.posre = f.read()

        with open(PATHS_EXAMPLEDATA.MMPBSA_LIG_ITP, "rb") as f:
            self.lig_itp = f.read()

        with open(PATHS_EXAMPLEDATA.MMPBSA_LIG_POSRE, "rb") as f:
            self.lig_posre = f.read()

        self.topol = GromacsState()
        top_path = os.path.dirname(PATHS_EXAMPLEDATA.GROMACS_1BVG_TOP)
        top_file = PATHS_EXAMPLEDATA.GROMACS_1BVG_TOP.split("/")[-1]
        self.topol.parse(top_path, top_file)
        self.topol.structures = {
            0: GenericData(_SGE.STD_STRUCTURE, file_data=self.structure)
        }
        self.topol.add_itp(
            os.path.join(MAIN_CONFIG["ICOLOS_TEST_DATA"], "gromacs/protein"),
            ["DMP:100.itp"],
        )
        self.topol.trajectories = {
            0: GenericData(file_name=_SGE.STD_XTC, file_data=self.xtc_file)
        }
        self.topol.tprs = {0: GenericData(_SGE.STD_TPR, self.tpr_file)}

    def test_protein_lig_single_traj(self):
        step_conf = {
            _SBE.STEPID: "test_gmmpbsa",
            _SBE.STEP_TYPE: "gmx_mmpbsa",
            _SBE.EXEC: {_SBE.EXEC_PREFIXEXECUTION: "module load gmx_MMPBSA "},
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
        step_mmpbsa.data.gmx_state = self.topol
        step_mmpbsa.data.generic.add_file(
            GenericData(file_name="DMP:100.itp", file_data=self.lig_itp)
        )
        step_mmpbsa.execute()
        out_path = os.path.join(self._test_dir, "ICOLOS_PROPS.dat")
        step_mmpbsa.data.gmx_state.write_props(self._test_dir)
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 200)

    def test_protein_lig_single_traj_custom_file(self):

        step_conf = {
            _SBE.STEPID: "test_gmmpbsa",
            _SBE.STEP_TYPE: "gmx_mmpbsa",
            _SBE.EXEC: {_SBE.EXEC_PREFIXEXECUTION: "module load gmx_MMPBSA"},
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                },
                _SBE.SETTINGS_ADDITIONAL: {
                    _SGE.FORCEFIELD: MAIN_CONFIG["FORCEFIELD"],
                    _SGE.COUPLING_GROUPS: "Protein Other",
                    _SGE.INPUT_FILE: PATHS_EXAMPLEDATA.MMPBSA_CUSTOM_INPUT,
                },
            },
        }
        step_mmpbsa = StepGMXmmpbsa(**step_conf)
        step_mmpbsa.data.gmx_state = self.topol

        step_mmpbsa.data.generic.add_file(
            GenericData(file_name="DMP:100.itp", file_data=self.lig_itp)
        )
        step_mmpbsa.execute()
        out_path = os.path.join(self._test_dir, "ICOLOS_PROPS.dat")
        step_mmpbsa.data.gmx_state.write_props(self._test_dir)
        stat_inf = os.stat(out_path)

        self.assertGreater(stat_inf.st_size, 180)

    def test_protein_lig_single_traj_MPI(self):

        step_conf = {
            _SBE.STEPID: "test_gmmpbsa",
            _SBE.STEP_TYPE: "gmx_mmpbsa",
            _SBE.EXEC: {_SBE.EXEC_PREFIXEXECUTION: "module load gmx_MMPBSA"},
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: ["MPI"],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                },
                _SBE.SETTINGS_ADDITIONAL: {
                    _SGE.FORCEFIELD: MAIN_CONFIG["FORCEFIELD"],
                    _SGE.COUPLING_GROUPS: "Protein Other",
                    _SGE.INPUT_FILE: PATHS_EXAMPLEDATA.MMPBSA_CUSTOM_INPUT,
                    _SGE.THREADS: 4,
                },
            },
        }
        step_mmpbsa = StepGMXmmpbsa(**step_conf)
        step_mmpbsa.data.gmx_state = self.topol

        step_mmpbsa.data.generic.add_file(
            GenericData(file_name="DMP:100.itp", file_data=self.lig_itp)
        )
        step_mmpbsa.execute()
        out_path = os.path.join(self._test_dir, "ICOLOS_PROPS.dat")
        step_mmpbsa.data.gmx_state.write_props(self._test_dir)
        stat_inf = os.stat(out_path)

        self.assertGreater(stat_inf.st_size, 180)

    def test_protein_lig_multi_traj_MPI(self):
        step_conf = {
            _SBE.STEPID: "test_gmmpbsa",
            _SBE.STEP_TYPE: "gmx_mmpbsa",
            _SBE.EXEC: {_SBE.EXEC_PREFIXEXECUTION: "module load gmx_MMPBSA "},
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: ["MPI"],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                },
                _SBE.SETTINGS_ADDITIONAL: {
                    _SGE.FORCEFIELD: MAIN_CONFIG["FORCEFIELD"],
                    _SGE.COUPLING_GROUPS: "Protein Other",
                    _SGE.THREADS: 4,
                },
            },
        }
        topol = self.topol
        topol.structures[1] = GenericData(_SGE.STD_STRUCTURE, file_data=self.structure)
        topol.structures[2] = GenericData(_SGE.STD_STRUCTURE, file_data=self.structure)
        topol.structures[3] = GenericData(_SGE.STD_STRUCTURE, file_data=self.structure)
        topol.trajectories[1] = GenericData(
            file_name=_SGE.STD_XTC, file_data=self.xtc_file
        )
        topol.trajectories[2] = GenericData(
            file_name=_SGE.STD_XTC, file_data=self.xtc_file
        )
        topol.trajectories[3] = GenericData(
            file_name=_SGE.STD_XTC, file_data=self.xtc_file
        )
        topol.tprs[1] = GenericData(_SGE.STD_TPR, self.tpr_file)
        topol.tprs[2] = GenericData(_SGE.STD_TPR, self.tpr_file)
        topol.tprs[3] = GenericData(_SGE.STD_TPR, self.tpr_file)

        step_mmpbsa = StepGMXmmpbsa(**step_conf)
        step_mmpbsa.data.generic.add_file(
            GenericData(file_name="DMP:100.itp", file_data=self.lig_itp)
        )
        step_mmpbsa.data.gmx_state = topol
        step_mmpbsa.execute()
        out_path = os.path.join(self._test_dir, "ICOLOS_PROPS.dat")
        step_mmpbsa.data.gmx_state.write_props(self._test_dir)
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 180)
        self.assertEqual(
            np.mean(step_mmpbsa.get_topol().properties[_SGE.MMGBSA_DG]), -32.1691
        )
