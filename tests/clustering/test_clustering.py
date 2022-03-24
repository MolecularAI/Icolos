import unittest

from icolos.core.containers.compound import Compound, Enumeration
from icolos.core.workflow_steps.calculation.clustering import StepClustering

from icolos.utils.enums.step_enums import StepBaseEnum, StepClusteringEnum

from tests.tests_paths import PATHS_EXAMPLEDATA, get_mol_as_Conformer

_SBE = StepBaseEnum
_SC = StepClusteringEnum()


class Test_Clustering(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass

    def setUp(self):
        pass

    @classmethod
    def tearDownClass(cls):
        pass

    def test_Clustering(self):
        step_conf = {
            _SBE.STEPID: "01_clustering",
            _SBE.STEP_TYPE: _SBE.STEP_CLUSTERING,
            _SBE.SETTINGS: {
                _SBE.SETTINGS_ARGUMENTS: {
                    _SBE.SETTINGS_ARGUMENTS_FLAGS: [],
                    _SBE.SETTINGS_ARGUMENTS_PARAMETERS: {
                        _SC.N_CLUSTERS: 2,
                        _SC.MAX_ITER: 300,
                    },
                },
                _SBE.SETTINGS_ADDITIONAL: {
                    _SC.TOP_N_PER_SOLVENT: 3,
                    _SC.FEATURES: ["area", "dipole", "HB_acc"],
                    _SC.FREE_ENERGY_SOLVENT_TAGS: [
                        "G_h2o",
                        "G_meoh",
                        "G_octanol",
                        # "G_dmso",
                        # "G_cychex",
                        # "G_acn",
                        # "G_thf",
                    ],
                },
            },
        }

        cl_step = StepClustering(**step_conf)
        cl_step.get_compounds().append(Compound(compound_number=1))
        cl_step.get_compounds()[0].add_enumeration(Enumeration(), auto_update=True)
        conformers = get_mol_as_Conformer(PATHS_EXAMPLEDATA.CLUSTERING_11CONFS)
        cl_step.data.compounds[0][0].add_conformers(conformers, auto_update=True)

        # 11 conformers are put in, but due to clustering only 10 should come out; note, that if only one solvent was
        # selected, only 9 would be outputted (this is because 2 of the clusters have 4 members and TOP_N_PER_SOLVENT
        # is set to 3)
        self.assertEqual(len(cl_step.get_compounds()[0][0].get_conformers()), 11)
        cl_step.execute()
        self.assertEqual(len(cl_step.get_compounds()[0][0].get_conformers()), 10)

        # make sure it is the 10th element (index 9) that has been removed
        self.assertListEqual(
            [
                list(
                    cl_step.get_compounds()[0][0]
                    .get_conformers()[i]
                    .get_molecule()
                    .GetConformer(0)
                    .GetPositions()[0]
                )
                for i in range(10)
            ],
            [
                [5.3347, 12.9328, 24.6745],
                [7.3326, 14.1647, 24.7427],
                [4.6181, 14.1638, 24.8496],
                [5.8612, 12.253, 24.0682],
                [5.1326, 14.4188, 25.3412],
                [6.6244, 13.1353, 25.9197],
                [5.4787, 13.5311, 24.3464],
                [7.3929, 12.7985, 23.9242],
                [5.7113, 13.2876, 24.6473],
                [7.3285, 12.8361, 23.9048],
            ],
        )
