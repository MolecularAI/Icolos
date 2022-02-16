import unittest
from copy import deepcopy

from rdkit.Geometry.rdGeometry import Point3D

from icolos.core.composite_agents.workflow import WorkFlow
from icolos.core.containers.compound import Compound, Conformer, Enumeration
from tests.tests_paths import PATHS_1UYD
from icolos.core.containers.generic import GenericData
from icolos.core.workflow_steps.io.data_manipulation import StepDataManipulation
from icolos.core.workflow_steps.step import StepBase
from icolos.utils.enums.composite_agents_enums import WorkflowEnum
from rdkit.Chem import SDMolSupplier

from icolos.utils.enums.step_enums import (
    StepBaseEnum,
    StepDataManipulationEnum,
    StepFilterEnum,
)
from icolos.utils.general.files_paths import attach_root_path, empty_output_dir
import os
from tests.tests_paths import (
    PATHS_EXAMPLEDATA,
    export_unit_test_env_vars,
    get_ligands_as_compounds_with_conformers,
    get_mol_as_Conformer,
)

_SBE = StepBaseEnum
_SDM = StepDataManipulationEnum()
_WE = WorkflowEnum()
_SFE = StepFilterEnum()


class Test_DataManipulation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/data_manip")
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

        export_unit_test_env_vars()

    def _get_step_filter_compounds(self):
        # produce the compounds object for testing
        # return 5 dummy compounds with 1 enumeration and 5 conformers per enumeration
        mols = SDMolSupplier(attach_root_path(PATHS_1UYD.LIGANDS))
        mol = mols[0]
        compounds = []
        for i in range(5):
            compound = Compound(name=str(i), compound_number=i)
            enum = Enumeration()
            for i in range(5):
                conf = Conformer(conformer=mol, conformer_id=i, enumeration_object=enum)
                enum.add_conformer(conformer=conf)
            compound.add_enumeration(enum)
            compounds.append(compound)
        return compounds

    def setUp(self):
        self._compounds = self._get_step_filter_compounds()
        # dG score gets gradually worse, prime gets gradually worse during the conformer walk
        dG_value = -13
        prime_value = -2900
        for compound in self._compounds:
            for enum in compound.get_enumerations():
                for conformer in enum.get_conformers():
                    conformer.get_molecule().SetProp("dG", str(dG_value))
                    conformer.get_molecule().SetProp(
                        "r_psp_MMGBSA_dG_Bind", str(prime_value)
                    )
                    dG_value += 0.2
                    prime_value -= 10
        # remove files from previous runs
        empty_output_dir(self._test_dir)

        with open(PATHS_1UYD.APO_MAE, "r") as f:
            data = f.read()
        self.complex_conformers = get_ligands_as_compounds_with_conformers(
            PATHS_1UYD.LIG4_POSES
        )
        self.mae_file = GenericData(file_name="structure.mae", file_data=data)
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

        # workflow (necessary to pass on data)
        conf = {
            _WE.HEADER: {_WE.ID: "test_workflow", _WE.ENVIRONMENT: {}},
            _WE.STEPS: [],
        }
        self.workflow = WorkFlow(**conf)

        # dummy step
        step_conf = {
            _SBE.STEPID: "01_dummy",
            _SBE.STEP_TYPE: "dummy",
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                },
                _SBE.SETTINGS_ADDITIONAL: {},
            },
        }
        dummy_step = StepBase(**step_conf)
        dummy_step.get_compounds().append(self.comp2)
        dummy_step.set_workflow_object(self.workflow)
        self.workflow.add_step(dummy_step)

    @classmethod
    def tearDownClass(cls):
        pass

    def test_matching(self):
        step_conf = {
            _SBE.STEPID: "01_data_manip",
            _SBE.STEP_TYPE: _SBE.STEP_DATA_MANIPULATION,
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                },
                _SBE.SETTINGS_ADDITIONAL: {
                    _SDM.ACTION: _SDM.ACTION_ATTACH_CONFORMERS_AS_EXTRA
                },
            },
        }
        manip_step = StepDataManipulation(**step_conf)
        manip_step.get_compounds().append(self.comp1)
        manip_step.set_workflow_object(self.workflow)
        self.workflow.add_step(manip_step)

        manip_step.settings.additional[_SDM.MATCH_SOURCE] = "01_dummy"
        manip_step.execute()

        self.assertEqual(
            manip_step.get_compounds()[0][0][2]
            .get_extra_data()[_SDM.KEY_MATCHED][0]
            .get_index_string(),
            "1:0:2",
        )

    def test_no_action(self):
        step_conf = {
            _SBE.STEPID: "01_data_manip",
            _SBE.STEP_TYPE: _SBE.STEP_DATA_MANIPULATION,
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                },
                _SBE.SETTINGS_ADDITIONAL: {_SDM.ACTION: _SDM.ACTION_NO_ACTION},
            },
        }
        manip_step = StepDataManipulation(**step_conf)
        manip_step.get_compounds().append(self.comp1)
        manip_step.set_workflow_object(self.workflow)
        self.workflow.add_step(manip_step)

        manip_step.settings.additional[_SDM.MATCH_SOURCE] = "01_dummy"
        manip_step.execute()

        self.assertEqual(len(manip_step.get_compounds()[0][0]), 11)

    def test_convert_mae2pdb(self):
        step_conf = {
            _SBE.STEPID: "01_data_manip",
            _SBE.STEP_TYPE: _SBE.STEP_DATA_MANIPULATION,
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                },
                _SBE.SETTINGS_ADDITIONAL: {_SDM.ACTION: _SDM.CONVERT_MAE_TO_PDB},
            },
        }
        manip_step = StepDataManipulation(**step_conf)
        manip_step.set_workflow_object(self.workflow)
        manip_step.data.generic.add_file(self.mae_file)
        self.workflow.add_step(manip_step)

        manip_step.execute()
        out_path = os.path.join(self._test_dir, "structure.pdb")
        manip_step.write_generic_by_extension(self._test_dir, "pdb")
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 201300)

    def test_get_complexes(self):
        step_conf = {
            _SBE.STEPID: "01_data_manip",
            _SBE.STEP_TYPE: _SBE.STEP_DATA_MANIPULATION,
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {},
                },
                _SBE.SETTINGS_ADDITIONAL: {
                    _SDM.ACTION: _SDM.ASSEMBLE_COMPLEXES,
                    _SDM.RECEPTOR: PATHS_1UYD.PDB_PATH,
                },
            },
        }
        manip_step = StepDataManipulation(**step_conf)
        manip_step.data.compounds = self.complex_conformers
        manip_step.set_workflow_object(self.workflow)
        self.workflow.add_step(manip_step)

        manip_step.execute()
        out_path = os.path.join(self._test_dir, "0:0:0.pdb")
        manip_step.write_generic_by_extension(self._test_dir, "pdb")
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 6600)

    def test_filtering(self):
        step_conf = {
            _SBE.STEPID: "01_filtering",
            _SBE.STEP_TYPE: _SBE.STEP_DATA_MANIPULATION,
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ADDITIONAL: {
                    _SDM.ACTION: _SDM.FILTER,
                    _SFE.FILTER_LEVEL: "enumerations",
                    _SFE.CRITERIA: "dG",
                    _SFE.RETURN_N: 3,
                    _SFE.HIGHEST_IS_BEST: False,
                }
            },
        }

        step_filter = StepDataManipulation(**step_conf)
        step_filter.data.compounds = self._compounds

        step_filter.execute()
        dG_max = (
            step_filter.data.compounds[0]
            .get_enumerations()[0]
            .get_conformers()[0]
            .get_molecule()
            .GetProp("dG")
        )
        step_filter.write_conformers(
            path=os.path.join(self._test_dir, "filtered_confs.sdf")
        )
        out_path = os.path.join(self._test_dir, "filtered_confs.sdf")
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 39708)
        self.assertEqual(int(dG_max), -13)

    def test_combined_filtering(self):
        # filter based on a sum of two criteria attached to each conformer
        step_conf = {
            _SBE.STEPID: "01_filtering",
            _SBE.STEP_TYPE: _SBE.STEP_DATA_MANIPULATION,
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ADDITIONAL: {
                    _SDM.ACTION: _SDM.FILTER,
                    _SFE.FILTER_LEVEL: "enumerations",
                    _SFE.CRITERIA: ["dG", "r_psp_MMGBSA_dG_Bind"],
                    _SFE.RETURN_N: 3,
                    _SFE.HIGHEST_IS_BEST: False,
                    _SFE.AGGREGATION: "sum",
                }
            },
        }

        step_filter = StepDataManipulation(**step_conf)
        step_filter.data.compounds = self._compounds
        step_filter.execute()

        dG_bind_max = (
            step_filter.data.compounds[0]
            .get_enumerations()[0]
            .get_conformers()[0]
            .get_molecule()
            .GetProp("r_psp_MMGBSA_dG_Bind")
        )
        # check we can get single values pack properly
        self.assertEqual(int(dG_bind_max), -2900)
        self.assertEqual(len(step_filter.data.compounds), 5)
        self.assertEqual(len(step_filter.data.compounds[0][0].get_conformers()), 3)

    def test_combined_filtering_compound_level(self):
        # filter at the compound level
        step_conf = {
            _SBE.STEPID: "01_filtering",
            _SBE.STEP_TYPE: _SBE.STEP_DATA_MANIPULATION,
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ADDITIONAL: {
                    _SDM.ACTION: _SDM.FILTER,
                    _SFE.FILTER_LEVEL: "compounds",
                    _SFE.CRITERIA: ["dG", "r_psp_MMGBSA_dG_Bind"],
                    _SFE.RETURN_N: 3,
                    _SFE.HIGHEST_IS_BEST: False,
                    _SFE.AGGREGATION: "sum",
                }
            },
        }

        step_filter = StepDataManipulation(**step_conf)
        step_filter.data.compounds = self._compounds
        step_filter.execute()

        dG_bind_max = (
            step_filter.data.compounds[0]
            .get_enumerations()[0]
            .get_conformers()[0]
            .get_molecule()
            .GetProp("r_psp_MMGBSA_dG_Bind")
        )
        # check we can get single values pack properly
        self.assertEqual(int(dG_bind_max), -2900)
