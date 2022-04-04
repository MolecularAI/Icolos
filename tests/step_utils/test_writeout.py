import os
import unittest
from icolos.core.composite_agents.workflow import WorkFlow, WorkflowData
from icolos.core.containers.generic import GenericContainer, GenericData
from icolos.core.containers.compound import Compound, Enumeration
from icolos.core.step_utils.input_preparator import StepData
from icolos.core.step_utils.step_writeout import WriteOutHandler
from icolos.utils.enums.step_enums import StepBaseEnum
from icolos.utils.general.files_paths import attach_root_path, remove_folder
from tests.tests_paths import PATHS_1UYD, PATHS_EXAMPLEDATA, get_mol_as_Conformer
import shutil

_SBE = StepBaseEnum


class Test_WriteOut(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._test_dir = attach_root_path("tests/junk/write-out")
        remove_folder(cls._test_dir)
        if not os.path.isdir(cls._test_dir):
            os.makedirs(cls._test_dir)

    def setUp(self):
        comp = Compound(compound_number=1, name="paracetamol")
        enum_mol = get_mol_as_Conformer(PATHS_EXAMPLEDATA.PARACETAMOL_PATH)[
            0
        ].get_molecule()
        comp.add_enumeration(Enumeration(molecule=enum_mol), auto_update=True)
        conformers = get_mol_as_Conformer(PATHS_EXAMPLEDATA.CLUSTERING_11CONFS)
        comp[0].add_conformers(conformers, auto_update=True)
        self.compound = comp

        comp2 = Compound(compound_number=2)
        comp2.add_enumeration(Enumeration(molecule=enum_mol), auto_update=True)
        conformers = get_mol_as_Conformer(PATHS_EXAMPLEDATA.CLUSTERING_11CONFS)
        comp2[0].add_conformers(conformers, auto_update=True)
        self.compound2 = comp2

    @classmethod
    def tearDownClass(cls):
        pass

    def test_conformer_writeout_merged(self):
        conf = {
            _SBE.WRITEOUT_CONFIG: {
                _SBE.WRITEOUT_COMP: {
                    _SBE.WRITEOUT_COMP_CATEGORY: _SBE.WRITEOUT_COMP_CATEGORY_CONFORMERS
                },
                _SBE.WRITEOUT_DESTINATION: {
                    _SBE.WRITEOUT_DESTINATION_TYPE: _SBE.WRITEOUT_DESTINATION_TYPE_FILE,
                    _SBE.WRITEOUT_DESTINATION_FORMAT: _SBE.FORMAT_SDF,
                },
            }
        }
        writeout_handler = WriteOutHandler(**conf)
        writeout_handler.set_data(StepData(compounds=[self.compound, self.compound2]))

        writeout_handler.config.destination.resource = os.path.join(
            self._test_dir, "both_compounds.sdf"
        )
        writeout_handler.write()
        stat_inf = os.stat(os.path.join(self._test_dir, "both_compounds.sdf"))
        self.assertGreater(stat_inf.st_size, 32100)

    def test_conformer_writeout_split(self):
        conf = {
            _SBE.WRITEOUT_CONFIG: {
                _SBE.WRITEOUT_COMP: {
                    _SBE.WRITEOUT_COMP_CATEGORY: _SBE.WRITEOUT_COMP_CATEGORY_CONFORMERS
                },
                _SBE.WRITEOUT_DESTINATION: {
                    _SBE.WRITEOUT_DESTINATION_TYPE: _SBE.WRITEOUT_DESTINATION_TYPE_FILE,
                    _SBE.WRITEOUT_DESTINATION_FORMAT: _SBE.FORMAT_SDF,
                    _SBE.WRITEOUT_DESTINATION_MERGE: False,
                },
            }
        }
        writeout_handler = WriteOutHandler(**conf)
        writeout_handler.set_data(StepData(compounds=[self.compound, self.compound2]))

        writeout_handler.config.destination.resource = os.path.join(
            self._test_dir, "[compound_id]_split.sdf"
        )
        writeout_handler.write()
        stat_inf = os.stat(os.path.join(self._test_dir, "1_split.sdf"))
        self.assertGreater(stat_inf.st_size, 16200)
        stat_inf = os.stat(os.path.join(self._test_dir, "2_split.sdf"))
        self.assertGreater(stat_inf.st_size, 15800)

    def test_extradata_writeout(self):
        conf = {
            _SBE.WRITEOUT_CONFIG: {
                _SBE.WRITEOUT_COMP: {
                    _SBE.WRITEOUT_COMP_CATEGORY: _SBE.WRITEOUT_COMP_CATEGORY_EXTRADATA,
                    _SBE.WRITEOUT_COMP_KEY: "testdata",
                },
                _SBE.WRITEOUT_DESTINATION: {
                    _SBE.WRITEOUT_DESTINATION_TYPE: _SBE.WRITEOUT_DESTINATION_TYPE_FILE,
                    _SBE.WRITEOUT_DESTINATION_FORMAT: _SBE.FORMAT_TXT,
                },
            }
        }
        self.compound[0][0].add_extra_data("testdata", ["this\n", "is\n", "a\ntest"])
        self.compound[0][1].add_extra_data(
            "testdata", "YETANOTHERTEST\nthis\nis\na\ntest"
        )
        self.compound[0]._conformers = [self.compound[0][0], self.compound[0][1]]

        writeout_handler = WriteOutHandler(**conf)
        writeout_handler.set_data(StepData(compounds=[self.compound]))

        # generate two files
        writeout_handler.config.destination.resource = os.path.join(
            self._test_dir, "extra_writeout/[conformer_id].txt"
        )
        writeout_handler.write()
        stat_inf = os.stat(os.path.join(self._test_dir, "extra_writeout/0.txt"))
        self.assertEqual(stat_inf.st_size, 15)
        stat_inf = os.stat(os.path.join(self._test_dir, "extra_writeout/1.txt"))
        self.assertEqual(stat_inf.st_size, 29)

    def test_tabular_writeout(self):
        config = {
            _SBE.WRITEOUT_CONFIG: {
                _SBE.WRITEOUT_COMP: {
                    _SBE.WRITEOUT_COMP_CATEGORY: _SBE.WRITEOUT_COMP_CATEGORY_CONFORMERS
                },
                _SBE.WRITEOUT_DESTINATION: {
                    _SBE.WRITEOUT_DESTINATION_TYPE: _SBE.WRITEOUT_DESTINATION_TYPE_FILE,
                    _SBE.WRITEOUT_DESTINATION_FORMAT: _SBE.FORMAT_CSV,
                },
            }
        }
        for idx, conf in enumerate(self.compound[0].get_conformers()):
            conf.get_molecule().SetProp("Gsolv_whatever", str(idx))
        self.compound[0][3].get_molecule().SetProp("another_prop", "bbc")
        writeout_handler = WriteOutHandler(**config)
        writeout_handler.set_data(StepData(compounds=[self.compound]))

        # write-out without selecting any tags
        writeout_handler.config.destination.resource = os.path.join(
            self._test_dir, "tabular_notagsselected_[conformer_id].csv"
        )
        writeout_handler.write()
        stat_inf = os.stat(os.path.join(self._test_dir, "tabular_notagsselected_0.csv"))
        self.assertGreater(stat_inf.st_size, 250)

        # write-out with selecting tags
        writeout_handler.config.destination.resource = os.path.join(
            self._test_dir, "tabular_tagsselected_[conformer_id].csv"
        )
        writeout_handler.config.compounds.selected_tags = [
            "Gsolv_whatever",
            "Gsolv_dmso",
            "another_prop",
        ]
        writeout_handler.write()
        stat_inf = os.stat(os.path.join(self._test_dir, "tabular_tagsselected_0.csv"))
        self.assertGreater(stat_inf.st_size, 300)

    def test_tabular_writeout_aggregate(self):
        config = {
            _SBE.WRITEOUT_CONFIG: {
                _SBE.WRITEOUT_COMP: {
                    _SBE.WRITEOUT_COMP_CATEGORY: _SBE.WRITEOUT_COMP_CATEGORY_CONFORMERS
                },
                _SBE.WRITEOUT_DESTINATION: {
                    _SBE.WRITEOUT_DESTINATION_TYPE: _SBE.WRITEOUT_DESTINATION_TYPE_FILE,
                    _SBE.WRITEOUT_DESTINATION_FORMAT: _SBE.FORMAT_CSV,
                },
            }
        }
        for idx, conf in enumerate(self.compound[0].get_conformers()):
            conf.get_molecule().SetProp("Gsolv_whatever", str(idx))
        self.compound[0][3].get_molecule().SetProp("another_prop", "bbc")
        writeout_handler = WriteOutHandler(**config)
        writeout_handler.set_data(StepData(compounds=[self.compound]))

        # write-out without selecting tags and using compound-level aggregation
        writeout_handler.config.destination.resource = os.path.join(
            self._test_dir, "tabular_notagsselected_[conformer_id]_compagg.csv"
        )
        writeout_handler.config.compounds.aggregation.mode = (
            _SBE.WRITEOUT_COMP_AGGREGATION_MODE_BESTPERCOMPOUND
        )
        writeout_handler.config.compounds.selected_tags = ["area"]
        writeout_handler.config.compounds.aggregation.key = "area"
        writeout_handler.write()
        stat_inf = os.stat(
            os.path.join(self._test_dir, "tabular_notagsselected_8_compagg.csv")
        )
        self.assertEqual(stat_inf.st_size, 61)

        # write-out without selecting tags and using compound-level aggregation (reverse)
        writeout_handler.config.destination.resource = os.path.join(
            self._test_dir, "tabular_notagsselected_[conformer_id]_compagg.csv"
        )
        writeout_handler.config.compounds.aggregation.mode = (
            _SBE.WRITEOUT_COMP_AGGREGATION_MODE_BESTPERCOMPOUND
        )
        writeout_handler.config.compounds.selected_tags = ["area"]
        writeout_handler.config.compounds.aggregation.key = "area"
        writeout_handler.config.compounds.aggregation.highest_is_best = False
        writeout_handler.write()
        stat_inf = os.stat(
            os.path.join(self._test_dir, "tabular_notagsselected_0_compagg.csv")
        )
        self.assertEqual(stat_inf.st_size, 61)

    def test_reinvent_writeout_empty(self):
        config = {
            _SBE.WRITEOUT_CONFIG: {
                _SBE.WRITEOUT_COMP: {
                    _SBE.WRITEOUT_COMP_CATEGORY: _SBE.WRITEOUT_COMP_CATEGORY_CONFORMERS
                },
                _SBE.WRITEOUT_DESTINATION: {
                    _SBE.WRITEOUT_DESTINATION_TYPE: _SBE.WRITEOUT_DESTINATION_TYPE_REINVENT,
                    _SBE.WRITEOUT_DESTINATION_FORMAT: _SBE.FORMAT_JSON,
                    _SBE.WRITEOUT_DESTINATION_RESOURCE: os.path.join(
                        self._test_dir, "reinvent_empty.json"
                    ),
                },
            }
        }
        for idx, conf in enumerate(self.compound[0].get_conformers()):
            conf.get_molecule().SetProp("Gsolv_whatever", str(idx))
        self.compound[0].clear_conformers()
        writeout_handler = WriteOutHandler(**config)
        writeout_handler.set_data(StepData(compounds=[self.compound]))

        writeout_handler.config.compounds.selected_tags = [
            "conformer_energy",
            "G_octanol",
        ]

        # write-out to console (REINVENT style)
        writeout_handler.write()

        # write-out to file
        out_path = os.path.join(self._test_dir, "reinvent_empty.json")
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 300)

    def test_reinvent_writeout_merged(self):
        config = {
            _SBE.WRITEOUT_CONFIG: {
                _SBE.WRITEOUT_COMP: {
                    _SBE.WRITEOUT_COMP_CATEGORY: _SBE.WRITEOUT_COMP_CATEGORY_CONFORMERS
                },
                _SBE.WRITEOUT_DESTINATION: {
                    _SBE.WRITEOUT_DESTINATION_TYPE: _SBE.WRITEOUT_DESTINATION_TYPE_REINVENT,
                    _SBE.WRITEOUT_DESTINATION_FORMAT: _SBE.FORMAT_JSON,
                    _SBE.WRITEOUT_DESTINATION_RESOURCE: os.path.join(
                        self._test_dir, "reinvent.json"
                    ),
                },
            }
        }
        for idx, conf in enumerate(self.compound[0].get_conformers()):
            conf.get_molecule().SetProp("Gsolv_whatever", str(idx))
        self.compound[0][3].get_molecule().SetProp("another_prop", "bbc")
        writeout_handler = WriteOutHandler(**config)
        writeout_handler.set_data(StepData(compounds=[self.compound]))

        writeout_handler.config.compounds.selected_tags = [
            "conformer_energy",
            "G_octanol",
        ]

        # write-out to console (REINVENT style)
        writeout_handler.write()

        # write-out to file
        out_path = os.path.join(self._test_dir, "reinvent.json")
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 300)

    def test_generic_writeout(self):
        conf = {
            _SBE.WRITEOUT_CONFIG: {
                _SBE.WRITEOUT_GENERIC: {_SBE.WRITEOUT_GENERIC_KEY: "txt"},
                _SBE.WRITEOUT_DESTINATION: {
                    _SBE.WRITEOUT_DESTINATION_TYPE: _SBE.WRITEOUT_DESTINATION_TYPE_FILE,
                    _SBE.WRITEOUT_DESTINATION_FORMAT: _SBE.FORMAT_TXT,
                },
            }
        }
        gc = GenericContainer()
        gc.add_file(
            GenericData(
                file_name="anothertest.txt",
                file_data="YETANOTHERTEST\nthis\nis\na\ntest",
            )
        )
        writeout_handler = WriteOutHandler(**conf)
        writeout_handler.set_data(StepData(generic=gc))

        # generate two files
        out_path = os.path.join(self._test_dir, "anothertest.txt")
        writeout_handler.config.destination.resource = out_path
        writeout_handler.write()
        out_path = os.path.join(self._test_dir, "anothertest_0.txt")
        stat_inf = os.stat(out_path)
        self.assertEqual(stat_inf.st_size, 29)

    def test_generic_writeout_path(self):
        conf = {
            _SBE.WRITEOUT_CONFIG: {
                _SBE.WRITEOUT_GENERIC: {_SBE.WRITEOUT_GENERIC_KEY: "xtc"},
                _SBE.WRITEOUT_DESTINATION: {
                    _SBE.WRITEOUT_DESTINATION_TYPE: _SBE.WRITEOUT_DESTINATION_TYPE_FILE,
                    _SBE.WRITEOUT_DESTINATION_FORMAT: _SBE.FORMAT_TXT,
                },
            }
        }

        writeout_handler = WriteOutHandler(**conf)
        # simulate the data being a path to a large file on disk
        gc = GenericContainer()
        gc.add_file(
            GenericData(
                file_name="md_0_1.xtc", file_data=PATHS_EXAMPLEDATA.GROMACS_1BVG_XTC
            )
        )
        writeout_handler.set_data(StepData(generic=gc))
        out_path = os.path.join(self._test_dir, "md_0_1.xtc")
        writeout_handler.config.destination.resource = out_path
        out_path = os.path.join(self._test_dir, "md_0_1_0.xtc")
        writeout_handler.write()

        # reset the files since by default it gets removed from the source location
        if not os.path.isfile(PATHS_EXAMPLEDATA.GROMACS_1BVG_XTC):
            shutil.copyfile(out_path, PATHS_EXAMPLEDATA.GROMACS_1BVG_XTC)
        stat_inf = os.stat(out_path)
        self.assertGreater(stat_inf.st_size, 16415800)

    def test_gmx_state_writeout(self):
        conf = {
            _SBE.WRITEOUT_CONFIG: {
                _SBE.WRITEOUT_GMX: {_SBE.WRITEOUT_GENERIC_KEY: "xtc, top, gro"},
                _SBE.WRITEOUT_DESTINATION: {
                    _SBE.WRITEOUT_DESTINATION_RESOURCE: self._test_dir
                },
            }
        }
        writeout_handler = WriteOutHandler(**conf)
        writeout_handler.data = StepData()

        writeout_handler.data.gmx_state.trajectories[0] = GenericData(
            file_name="traj.xtc", file_data=PATHS_EXAMPLEDATA.GROMACS_1BVG_XTC
        )
        writeout_handler.data.gmx_state.trajectories[1] = GenericData(
            file_name="traj.xtc", file_data=PATHS_EXAMPLEDATA.GROMACS_1BVG_XTC
        )
        writeout_handler.data.gmx_state.structures[0] = GenericData(
            file_name="confout.gro",
            file_data=PATHS_EXAMPLEDATA.GROMACS_HOLO_STRUCTURE_GRO,
        )
        writeout_handler.data.gmx_state.structures[1] = GenericData(
            file_name="confout.gro",
            file_data=PATHS_EXAMPLEDATA.GROMACS_HOLO_STRUCTURE_GRO,
        )
        writeout_handler.data.gmx_state.set_topol(
            path="", file=PATHS_EXAMPLEDATA.GROMACS_1BVG_TOP
        )

        writeout_handler.write()

        for file in [
            "confout_0.gro",
            "confout_1.gro",
            "topol.top",
            "traj_0.xtc",
            "traj_1.xtc",
        ]:
            assert os.path.isfile(os.path.join(self._test_dir, file))
