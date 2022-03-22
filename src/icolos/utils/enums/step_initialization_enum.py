from icolos.core.step_dispatch.dispatcher import StepDispatcher
from icolos.core.workflow_steps.autodockvina.docking import StepAutoDockVina
from icolos.core.workflow_steps.autodockvina.target_preparation import (
    StepAutoDockVinaTargetPreparation,
)
from icolos.core.workflow_steps.calculation.jazzy import StepJazzy
from icolos.core.workflow_steps.calculation.kallisto import StepKallisto
from icolos.core.workflow_steps.ccdc.docking import StepGold
from icolos.core.workflow_steps.calculation.electrostatics.esp_sim import StepEspSim
from icolos.core.workflow_steps.calculation.feature_counter import StepFeatureCounter
from icolos.core.workflow_steps.gromacs.do_dssp import StepGMXDoDSSP
from icolos.core.workflow_steps.gromacs.mmpbsa import StepGMXmmpbsa
from icolos.core.workflow_steps.fpocket.mdpocket import StepMDpocket
from icolos.core.workflow_steps.gromacs.trajcat import StepGMXTrjcat
from icolos.core.workflow_steps.io.data_manipulation import StepDataManipulation
from icolos.core.workflow_steps.schrodinger.fep_analysis import StepFepPlusAnalysis
from icolos.core.workflow_steps.structure_prediction.pdb_fixer import StepPdbFixer
from icolos.core.workflow_steps.gromacs import *
from icolos.core.workflow_steps.calculation.boltzmann_weighting import (
    StepBoltzmannWeighting,
)
from icolos.core.workflow_steps.calculation.rmsd import StepRMSD
from icolos.core.workflow_steps.schrodinger import *
from icolos.core.workflow_steps.calculation.cosmo import StepCosmo
from icolos.core.workflow_steps.calculation.turbomole import StepTurbomole
from icolos.core.workflow_steps.confgen.crest import StepCREST
from icolos.core.workflow_steps.pmx import *
from icolos.core.workflow_steps.confgen.omega import StepOmega
from icolos.core.workflow_steps.confgen.xtb import StepXTB
from icolos.core.workflow_steps.io.embedder import StepEmbedding
from icolos.core.workflow_steps.io.initialize_compound import StepInitializeCompound
from icolos.core.workflow_steps.prediction.predictor import StepPredictor
from icolos.core.workflow_steps.prediction.model_building import StepModelBuilder
from icolos.core.workflow_steps.calculation.clustering import StepClustering
from icolos.core.workflow_steps.calculation.rms_filter import StepRMSFilter
from icolos.core.workflow_steps.calculation.panther import StepPanther
from icolos.core.workflow_steps.calculation.shaep import StepShaep
from icolos.core.workflow_steps.structure_prediction.peptide_embedder import (
    StepPeptideEmbedder,
)
from icolos.core.workflow_steps.structure_prediction.dssp import StepDSSP
from icolos.utils.enums.step_enums import StepBaseEnum


_SBE = StepBaseEnum


class StepInitializationEnum:

    STEP_INIT_DICT = {
        _SBE.STEP_CREST: StepCREST,
        _SBE.STEP_OMEGA: StepOmega,
        _SBE.STEP_XTB: StepXTB,
        _SBE.STEP_MACROMODEL: StepMacromodel,
        _SBE.STEP_TURBOMOLE: StepTurbomole,
        _SBE.STEP_COSMO: StepCosmo,
        _SBE.STEP_INITIALIZATION: StepInitializeCompound,
        _SBE.STEP_EMBEDDING: StepEmbedding,
        _SBE.STEP_PREDICTION: StepPredictor,
        _SBE.STEP_MODEL_BUILDING: StepModelBuilder,
        _SBE.STEP_BOLTZMANN_WEIGHTING: StepBoltzmannWeighting,
        _SBE.STEP_PRIME: StepPrime,
        _SBE.STEP_DESMOND: StepDesmondExec,
        _SBE.STEP_DESMOND_SETUP: StepDesmondSetup,
        _SBE.STEP_CLUSTERING: StepClustering,
        _SBE.STEP_RMSFILTER: StepRMSFilter,
        _SBE.STEP_PANTHER: StepPanther,
        _SBE.STEP_KALLISTO: StepKallisto,
        _SBE.STEP_JAZZY: StepJazzy,
        _SBE.STEP_SHAEP: StepShaep,
        _SBE.STEP_PDB2GMX: StepGMXPdb2gmx,
        _SBE.STEP_EDITCONF: StepGMXEditConf,
        _SBE.STEP_SOLVATE: StepGMXSolvate,
        _SBE.STEP_GENION: StepGMXGenion,
        _SBE.STEP_GROMPP: StepGMXGrompp,
        _SBE.STEP_MDRUN: StepGMXMDrun,
        _SBE.STEP_TRJCONV: StepGMXTrjconv,
        _SBE.STEP_TRJCAT: StepGMXTrjcat,
        _SBE.STEP_CLUSTER: StepGMXCluster,
        _SBE.STEP_DO_DSSP: StepGMXDoDSSP,
        _SBE.STEP_GMX_RMSD: StepGMXrmsd,
        _SBE.STEP_LIGPREP: StepLigprep,
        _SBE.STEP_GLIDE: StepGlide,
        _SBE.STEP_FEP_PLUS_SETUP: StepFepPlusSetup,
        _SBE.STEP_FEP_PLUS_EXEC: StepFepPlusExec,
        _SBE.STEP_FEP_PLUS_ANALYSIS: StepFepPlusAnalysis,
        _SBE.STEP_PREPWIZARD: StepPrepwizard,
        _SBE.STEP_MDPOCKET: StepMDpocket,
        _SBE.STEP_PEPTIDE_EMBEDDER: StepPeptideEmbedder,
        _SBE.STEP_PDB_FIXER: StepPdbFixer,
        _SBE.STEP_GMX_MMPBSA: StepGMXmmpbsa,
        _SBE.STEP_TS_CLUSTER: StepClusterTS,
        _SBE.STEP_DSSP: StepDSSP,
        _SBE.STEP_RMSD: StepRMSD,
        _SBE.STEP_DATA_MANIPULATION: StepDataManipulation,
        _SBE.STEP_PMX_ASSEMBLE_SYSTEMS: StepPMXAssembleSystems,
        _SBE.STEP_PMX_ATOMMAPPING: StepPMXatomMapping,
        _SBE.STEP_PMX_ABFE_SETUP: StepPMXabfe,
        _SBE.STEP_PMX_DOUBLEBOX: StepPMXdoublebox,
        _SBE.STEP_PMX_LIGANDHYBRID: StepPMXligandHybrid,
        _SBE.STEP_PMX_BOX_WATER_IONS: StepPMXBoxWaterIons,
        _SBE.STEP_PMX_SETUP: StepPMXSetup,
        _SBE.STEP_PMX_PREPARE_TRANSITIONS: StepPMXPrepareTransitions,
        _SBE.STEP_PMX_PREPARE_SIMULATIONS: StepPMXPrepareSimulations,
        _SBE.STEP_PMX_RUN_ANALYSIS: StepPMXRunAnalysis,
        _SBE.STEP_PMX_MUTATE: StepPMXmutate,
        _SBE.STEP_PMX_GENTOP: StepPMXgentop,
        _SBE.STEP_PMX_GENLIB: StepPMXgenlib,
        _SBE.STEP_FEATURE_COUNTER: StepFeatureCounter,
        _SBE.STEP_AUTODOCKVINA_DOCKING: StepAutoDockVina,
        _SBE.STEP_AUTODOCKVINA_TARGET_PREPARATION: StepAutoDockVinaTargetPreparation,
        _SBE.STEP_GOLD_DOCKING: StepGold,
        _SBE.STEP_PMX_RUN_SIMULATIONS: StepPMXRunSimulations,
        _SBE.STEP_DISPATCHER: StepDispatcher,
        _SBE.STEP_ESP_SIM: StepEspSim,
    }
