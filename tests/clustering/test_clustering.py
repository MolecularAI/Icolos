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
                        _SC.N_CLUSTERS: 3,
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
                        "G_dmso",
                        "G_cychex",
                        "G_acn",
                        "G_thf",
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
                [0.8838, 0.6808, -0.1373],
                [-4.2269, -0.441, 0.2359],
                [-4.1693, -0.532, -0.0567],
                [-4.2326, -0.4701, 0.3534],
                [-4.201, -0.5446, 0.131],
                [-4.2198, -0.4705, 0.1656],
                [-4.2318, -0.444, 0.2474],
                [-4.2316, -0.14, 0.0848],
                [-4.1953, -0.1989, -0.1017],
                [-4.1528, -0.0208, 0.0932],
            ],
        )
