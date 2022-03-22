from enum import Enum


class StepBaseEnum(str, Enum):
    # general settings
    STEPID = "step_id"

    # different step types
    STEP_TYPE = "type"
    STEP_CREST = "CREST"
    STEP_OMEGA = "OMEGA"
    STEP_XTB = "XTB"
    STEP_MACROMODEL = "MACROMODEL"
    STEP_TURBOMOLE = "TURBOMOLE"
    STEP_COSMO = "COSMO"
    STEP_INITIALIZATION = "INITIALIZATION"
    STEP_EMBEDDING = "EMBEDDING"
    STEP_PREDICTION = "PREDICTION"
    STEP_MODEL_BUILDING = "MODEL_BUILDING"
    STEP_FEATURE_COUNTER = "FEATURE_COUNTER"
    STEP_BOLTZMANN_WEIGHTING = "BOLTZMANN_WEIGHTING"
    STEP_PKA_PREDICTION = "PKA_PREDICTION"
    STEP_PRIME = "PRIME"
    STEP_CLUSTERING = "CLUSTERING"
    STEP_RMSD = "RMSD"
    STEP_RMSFILTER = "RMSFILTER"
    STEP_DATA_MANIPULATION = "DATA_MANIPULATION"
    STEP_DESMOND = "DESMOND"
    STEP_DESMOND_SETUP = "DESMOND_SETUP"
    STEP_FILTER = "FILTER"
    STEP_PANTHER = "PANTHER"
    STEP_KALLISTO = "KALLISTO"
    STEP_JAZZY = "JAZZY"
    STEP_SHAEP = "SHAEP"
    STEP_PDB2GMX = "PDB2GMX"
    STEP_EDITCONF = "EDITCONF"
    STEP_SOLVATE = "SOLVATE"
    STEP_GENION = "GENION"
    STEP_GROMPP = "GROMPP"
    STEP_MDRUN = "MDRUN"
    STEP_TRJCONV = "TRJCONV"
    STEP_TRJCAT = "TRJCAT"
    STEP_GMX_RMSD = "GMX_RMSD"
    STEP_CLUSTER = "CLUSTER"
    STEP_DO_DSSP = "DO_DSSP"
    STEP_LIGPREP = "LIGPREP"
    STEP_GLIDE = "GLIDE"
    STEP_AUTODOCKVINA_DOCKING = "VINA_DOCKING"
    STEP_AUTODOCKVINA_TARGET_PREPARATION = "VINA_TARGET_PREPARATION"
    STEP_GOLD_DOCKING = "GOLD_DOCKING"
    STEP_FEP_PLUS_SETUP = "FEP_PLUS_SETUP"
    STEP_FEP_PLUS_EXEC = "FEP_PLUS_EXEC"
    STEP_FEP_PLUS_ANALYSIS = "FEP_PLUS_ANALYSIS"
    STEP_PREPWIZARD = "PREPWIZARD"
    STEP_MDPOCKET = "MDPOCKET"
    STEP_PDB_FIXER = "PDB_FIXER"
    STEP_PEPTIDE_EMBEDDER = "PEPTIDE_EMBEDDER"
    STEP_GMX_MMPBSA = "GMX_MMPBSA"

    # PMX SCRIPTS
    STEP_PMX_ABFE_SETUP = "PMX_ABFE_SETUP"
    STEP_PMX_ANALYSE = "PMX_ANALYSE"
    STEP_PMX_ATOMMAPPING = "PMX_ATOMMAPPING"
    STEP_PMX_DOUBLEBOX = "PMX_DOUBLEBOX"
    STEP_PMX_GENLIB = "PMX_GENLIB"
    STEP_PMX_GENTOP = "PMX_GENTOP"
    STEP_PMX_LIGANDHYBRID = "PMX_LIGANDHYBRID"
    STEP_PMX_MUTATE = "PMX_MUTATE"
    STEP_PMX_SETUP = "PMX_SETUP"
    STEP_PMX_PREPARE_SIMULATIONS = "PMX_PREPARE_SIMULATIONS"
    STEP_PMX_BOX_WATER_IONS = "PMX_BOX_WATER_IONS"
    STEP_PMX_PREPARE_TRANSITIONS = "PMX_PREPARE_TRANSITIONS"
    STEP_PMX_RUN_SIMULATIONS = "PMX_RUN_SIMULATIONS"
    STEP_PMX_ASSEMBLE_SYSTEMS = "PMX_ASSEMBLE_SYSTEMS"
    STEP_PMX_RUN_ANALYSIS = "PMX_RUN_ANALYSIS"

    STEP_DSSP = "DSSP"
    STEP_TS_CLUSTER = "TS_CLUSTER"
    STEP_ESP_SIM = "ESP_SIM"
    STEP_DISPATCHER = "DISPATCHER"
    STEP_ACTIVE_LEARNING = "ACTIVE_LEARNING"

    # flow control blocks
    STEP_ITERATOR = "ITERATOR"

    # execution
    EXEC = "execution"
    EXEC_PREFIXEXECUTION = "prefix_execution"
    EXEC_BINARYLOCATION = "binary_location"
    EXEC_PARALLELIZATION = "parallelization"
    EXEC_PARALLELIZATION_CORES = "jobs"
    EXEC_PARALLELIZATION_MAXLENSUBLIST = "max_length_sublists"
    EXEC_FAILUREPOLICY = "failure_policy"
    EXEC_FAILUREPOLICY_NTRIES = "n_tries"
    EXEC_PLATFORM = "platform"
    EXEC_RESOURCES = "resources"
    EXEC_RESOURCES_PARTITION = "partition"
    EXEC_RESOURCES_GRES = "gres"
    EXEC_RESOURCES_TASKS = "tasks"
    EXEC_RESOURCES_MODULES = "modules"
    EXEC_RESOURCES_MEM = "mem"
    EXEC_RESOURCES_CORES = "cores"
    EXEC_RESOURCES_OTHER_ARGS = "other_args"
    EXEC_RESOURCES_ADDITIONAL_LINES = "additional_lines"

    # settings
    SETTINGS = "settings"
    SETTINGS_ARGUMENTS = "arguments"
    SETTINGS_ARGUMENTS_FLAGS = "flags"
    SETTINGS_ARGUMENTS_PARAMETERS = "parameters"
    SETTINGS_ADDITIONAL = "additional"

    PIPE_INPUT = "pipe_input"

    # annotation: fixed strings
    ANNOTATION_TAG_DOCKING_SCORE = "docking_score"
    ANNOTATION_TAG_G_SCORE = "g_score"

    ANNOTATION_GRID_ID = "grid_id"
    ANNOTATION_GRID_PATH = "grid_path"
    ANNOTATION_GRID_FILENAME = "grid_filename"

    GRID_IDS = "grid_ids"  # enforces given list of IDs rather than indices

    # I/O and "hand-over"
    # ---------
    FORMAT_SDF = "SDF"
    FORMAT_CSV = "CSV"
    FORMAT_SMI = "SMI"
    FORMAT_MOL2 = "MOL2"
    FORMAT_XTB = "XTB"
    FORMAT_PDB = "PDB"
    FORMAT_PKL = "PKL"
    FORMAT_SMILES = "SMILES"
    FORMAT_PLAIN = "PLAIN"
    FORMAT_TXT = "TXT"
    FORMAT_JSON = "JSON"
    FORMAT_DTR = "DTR"
    FORMAT_CMS = "CMS"

    INPUT = "input"
    INPUT_FIELD = "field"
    INPUT_SOURCES = "sources"
    INPUT_COMPOUNDS = "compounds"
    INPUT_ENUMERATIONS = "enumerations"
    INPUT_EXTENSION = "extension"
    INPUT_SOURCE = "source"
    INPUT_GENERIC = "generic"
    INPUT_GMX = "gmx_state"
    INPUT_FORMAT = "format"
    INPUT_SOURCE_TYPE = "source_type"
    INPUT_SOURCE_TYPE_FILE = "file"
    INPUT_SOURCE_TYPE_DIR = "dir"
    INPUT_SOURCE_TYPE_PATH = "path"
    INPUT_SOURCE_TYPE_STEP = "step"
    INPUT_SOURCE_TYPE_STRING = "string"
    INPUT_SOURCE_TYPE_URL = "url"

    INPUT_ENFORCE_IDS = "enforce_ids"
    INPUT_ENFORCE_COMPOUND_IDS = "compound_ids"
    INPUT_ENFORCE_ENUMERATION_IDS = "enumeration_ids"

    INPUT_MERGE = "merge"
    INPUT_MERGE_COMPOUNDS = "compounds"
    INPUT_MERGE_COMPOUNDS_BY = "merge_compounds_by"
    INPUT_MERGE_ENUMERATIONS = "enumerations"
    INPUT_MERGE_ENUMERATIONS_BY = "merge_enumerations_by"
    INPUT_MERGE_BY_NAME = "name"
    INPUT_MERGE_BY_ID = "id"
    INPUT_MERGE_BY_SMILE = "smile"

    FILE_TYPE_PDB = "pdb"
    FILETYPE_TXT = "txt"
    FILE_SIZE_THRESHOLD = 2e9

    # CSV settings
    INPUT_CSV_DELIMITER = "delimiter"
    INPUT_CSV_COLUMNS = "columns"
    INPUT_CSV_SMILES_COLUMN = "smiles"
    INPUT_CSV_NAMES_COLUMN = "names"

    # write-out settings
    WRITEOUT = "writeout"
    WRITEOUT_CONFIG = "config"

    WRITEOUT_COMP = "compounds"
    WRITEOUT_COMP_CATEGORY = "category"
    WRITEOUT_COMP_CATEGORY_CONFORMERS = "conformers"
    WRITEOUT_COMP_CATEGORY_ENUMERATIONS = "enumerations"
    WRITEOUT_COMP_CATEGORY_EXTRADATA = "extra_data"
    WRITEOUT_COMP_KEY = "key"
    WRITEOUT_COMP_AGGREGATION = "aggregation"
    WRITEOUT_COMP_AGGREGATION_MODE = "mode"
    WRITEOUT_COMP_AGGREGATION_MODE_ALL = "all"
    WRITEOUT_COMP_AGGREGATION_MODE_BESTPERCOMPOUND = "best_per_compound"
    WRITEOUT_COMP_AGGREGATION_MODE_BESTPERENUMERATION = "best_per_enumeration"
    WRITEOUT_COMP_AGGREGATION_HIGHESTISBEST = "highest_is_best"
    WRITEOUT_COMP_AGGREGATION_KEY = "key"
    WRITEOUT_COMP_SELECTED_TAGS = "selected_tags"
    WRITEOUT_COMP_SELECTED_TAGS_KEY = "key"
    WRITEOUT_COMP_SELECTED_TAGS_HIGHESTISBEST = "highest_is_best"

    WRITEOUT_GENERIC = "generic"
    WRITEOUT_GENERIC_KEY = "key"

    WRITEOUT_GMX = "gmx_state"

    WRITEOUT_DESTINATION = "destination"
    WRITEOUT_DESTINATION_RESOURCE = "resource"
    WRITEOUT_DESTINATION_TYPE = "type"
    WRITEOUT_DESTINATION_TYPE_FILE = "file"
    WRITEOUT_DESTINATION_TYPE_REINVENT = "reinvent"
    WRITEOUT_DESTINATION_TYPE_STDOUT = "stdout"
    WRITEOUT_DESTINATION_TYPE_STDERR = "stderr"
    WRITEOUT_DESTINATION_TYPE_REST = "rest"
    WRITEOUT_DESTINATION_FORMAT = "format"
    WRITEOUT_DESTINATION_MERGE = "merge"
    WRITEOUT_DESTINATION_AUTOMATIC = "automatic"
    WRITEOUT_DESTINATION_BASE_NAME = "base_name"
    WRITEOUT_DESTINATION_DIR = "dir"
    WRITEOUT_DESTINATION_MODE = "mode"

    TOKEN_GUARD = "token_guard"

    # try to find the internal value and return
    # def __getattr__(self, name):
    #     if name in self:
    #         return name
    #     raise AttributeError

    # prohibit any attempt to set any values
    # def __setattr__(self, key, value):
    #     raise ValueError("No changes allowed.")


class IteratorEnum(str, Enum):
    N_ITERS = "n_iters"
    ALL = "all"
    SINGLE = "single"


class StepEmbeddingEnum:
    METHOD = "method"
    METHOD_RDKIT = "RDKIT"

    EMBED_AS = "embed_as"
    EMBED_AS_ENUMERATIONS = "enumerations"
    EMBED_AS_CONFORMERS = "conformers"

    RDKIT_PROTONATE = "protonate"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class StepClusteringEnum:
    N_CLUSTERS = "n_clusters"
    MAX_ITER = "max_iter"
    TOP_N_PER_SOLVENT = "top_n_per_solvent"
    FEATURES = "features"
    FREE_ENERGY_SOLVENT_TAGS = "free_energy_solvent_tags"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class StepRMSFilterEnum:
    THRESHOLD = "threshold"  # RMS threshold in Angstrom

    # order by this tag in picking the conformers
    ORDER_BY = "order_by"
    ORDER_ASCENDING = "ascending"

    METHOD = "method"  # calculation method
    METHOD_BEST = "best"  # RDkit's "GetBestRMS"
    METHOD_ALIGNMOL = "alignmol"  # RDkit's "AlignMol"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class StepRMSDEnum:

    METHOD = "method"  # calculation method
    METHOD_BEST = "best"  # RDkit's "GetBestRMS"
    METHOD_ALIGNMOL = "alignmol"  # RDkit's "AlignMol"

    RMSD_TAG = "rmsd"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class StepFeatureCounterEnum:

    FEATURE = "feature"
    LEVEL = "level"
    LEVEL_ENUMERATION = "enumeration"
    LEVEL_CONFORMER = "conformer"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class StepDataManipulationEnum:

    # specify actions that can be used
    ACTION = "action"
    ACTION_NO_ACTION = (
        # used to skip any calculation (e.g. to just pool input data)
        "no_action"
    )
    CONVERT_MAE_TO_PDB = "mae2pdb"
    ASSEMBLE_COMPLEXES = "assemble_complexes"
    ACTION_ATTACH_CONFORMERS_AS_EXTRA = "attach_conformers_as_extra"
    COLLECT_ITERATOR_RESULTS = "collect_iterator_results"
    FILTER = "filter"

    # --> For ACTION_ATTACH_CONFORMERS_AS_EXTRA
    # --- Match everything with the same <compound_id>:<enumeration_id>:<conformer_id> string
    MATCH_SOURCE = (
        "source"  # step from which the conformers are to be used for matching
    )
    KEY_MATCHED = "matched"  # extra-data key for matched data
    RECEPTOR = "receptor"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class StepBoltzmannWeightingEnum:
    PROPERTIES = "properties"
    PROPERTIES_INPUT = "input"
    PROPERTIES_OUTPUT = "output"

    WEIGHT = "weight"
    WEIGHT_INPUT = "input"
    WEIGHT_OUTPUT_PREFIX = "output_prefix"
    WEIGHT_PROPERTIES = "properties"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class StepPredictorEnum:
    MODEL_PATH = "model_path"
    FEATURES = "features"
    NAME_PREDICTED = "name_predicted"
    NAME_PREDICTED_DEFAULT = "pred_value"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class StepGoldEnum:

    CONFIGURATION = "configuration"
    GOLD_CONFIG_FILE = "gold_config_file"
    BLOCK_INDENT = "  "
    CONFIGURATION_START = "GOLD CONFIGURATION FILE"
    AUTOMATIC_SETTINGS = "AUTOMATIC SETTINGS"
    AUTOSCALE = "autoscale"
    POPULATION = "POPULATION"
    GENETIC_OPERATORS = "GENETIC OPERATORS"
    FLOOD_FILL = "FLOOD FILL"
    CAVITY_FILE = "cavity_file"
    DATA_FILES = "DATA FILES"
    LIGAND_DATA_FILE = "ligand_data_file"
    FLAGS = "FLAGS"
    TERMINATION = "TERMINATION"
    CONSTRAINTS = "CONSTRAINTS"
    COVALENT_BONDING = "COVALENT BONDING"
    SAVE_OPTIONS = "SAVE OPTIONS"
    FITNESS_FUNCTION_SETTINGS = "FITNESS FUNCTION SETTINGS"
    GOLD_FITFUNC_PATH = "gold_fitfunc_path"
    PROTEIN_DATA = "PROTEIN DATA"
    PROTEIN_DATAFILE = "protein_datafile"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class StepModelBuilderEnum:
    # configuration fields
    DATA = "data"
    DATA_INPUT_COLUMN = "input_column"
    DATA_RESPONSE_COLUMN = "response_column"
    DATA_TRAININGSET_FILE = "training_dataset_file"
    DATA_TESTSET_FILE = "test_dataset_file"

    # fixed tempfile names
    TMP_INPUT_CONFIG = "input_config.json"
    TMP_INPUT_DATA = "input_data.csv"
    TMP_OUTPUT_BEST_MODEL = "best_model.pkl"
    TMP_OUTPUT_BEST_PARAMETERS = "best_parameters.json"
    TMP_OUTPUT_PRODUCTION_MODEL = "production_model.pkl"

    # fields
    FIELD_KEY_PRODUCTION_MODEL = "production_model"
    FIELD_KEY_BEST_CONFIGURATION = "best_configuration"
    FIELD_KEY_INPUT_DATA = "input_data"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class TokenGuardEnum:
    TG = "token_guard"
    TG_PREFIX_EXECUTION = "prefix_execution"
    TG_BINARY_LOCATION = "binary_location"
    TG_TOKEN_POOLS = "token_pools"
    TG_WAIT_INTERVAL_SECONDS = "wait_interval_seconds"
    TG_WAIT_LIMIT_SECONDS = "wait_limit_seconds"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class StepPrepwizEnum:

    REMOVE_RES = "remove_res"
    LIGANDS = "ligands"
    COFACTOR_IDS = [
        "TDP",
        "FAD",
        "FMN",
        "NAD",
        "PNS",
        "COA",
        "PLP",
        "GSH",
        "BTN",
        "FFO",
        "B12",
        "ASC",
        "MQ7",
        "UQ1",
        "MGD",
        "H4B",
        "MDO",
        "SAM",
        "F43",
        "COM",
        "TP7",
        "HEA",
        "DPM",
        "PQQ",
        "TPQ",
        "TRQ",
        "LPA",
        "HEM",
    ]


class StepLigprepEnum:
    FILTER_FILE = "filter_file"

    # the SDF tag with <identifier>-# (where # is the number of the enumeration starting with '1')
    LIGPREP_VARIANTS = "s_lp_Variant"
    # number from 0 to 1 (sums up to 1 over all variants)
    LIGPREP_TAUTOMER_PROBABILITY = "r_lp_tautomer_probability"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class StepGlideEnum:
    # input specification parameters
    CONFIGURATION = "configuration"
    TIME_LIMIT_PER_TASK = "time_limit_per_task"
    MAESTRO_IN_FILE = "maestro_in_file"
    MAESTRO_IN_FILE_PATH = "path"

    # glide: fixed strings
    # the docking score (including "Epik" corrections")
    GLIDE_DOCKING_SCORE = "r_i_docking_score"
    # the "docking score" without "Epik" corrections
    GLIDE_GSCORE = "r_i_glide_gscore"
    # the index of the ligand in the input file (starting with '1')
    GLIDE_SOURCE_FILE_INDEX = "i_m_source_file_index"

    GLIDE_POSEVIEWER_FILE_KEY = "structures_pv.maegz"
    GLIDE_MAEGZ_DEFAULT_EXTENSION = "_pv.maegz"
    GLIDE_SDF_DEFAULT_EXTENSION = "_lib.sdfgz"
    GLIDE_LOG = ".log"
    GLIDE_SDF = ".sdf"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class StepMacromodelEnum:
    # COM file
    COM_FILE = "com_file"
    COM_FILE_PATH = "com_file.com"
    COM_FILE_DEFAULT = """ MMOD       0      1      0      0     0.0000     0.0000     0.0000     0.0000
 DEBG      55      0      0      0     0.0000     0.0000     0.0000     0.0000
 FFLD      16      1      0      0     1.0000     0.0000     0.0000     0.0000
 SOLV       3      1      0      0     0.0000     0.0000     0.0000     0.0000
 EXNB       0      0      0      0     0.0000     0.0000     0.0000     0.0000
 BDCO       0      0      0      0    89.4427 99999.0000     0.0000     0.0000
 READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000
 CRMS       0      0      0      0     0.0000     0.8000     0.0000     0.0000
 LMCS    1000      0      0      0     0.0000     0.0000     3.0000     6.0000
 NANT       0      0      0      0     0.0000     0.0000     0.0000     0.0000
 MCNV       1      5      0      0     0.0000     0.0000     0.0000     0.0000
 MCSS       2      0      0      0    27.0000     0.0000     0.0000     0.0000
 MCOP       1      0      0      0     0.5000     0.0000     0.0000     0.0000
 DEMX       0    833      0      0    27.0000    54.0000     0.0000     0.0000
 MSYM       0      0      0      0     0.0000     0.0000     0.0000     0.0000
 AUOP       0      0      0      0   100.0000     0.0000     0.0000     0.0000
 AUTO       0      2      1      1     0.0000     1.0000     0.0000     2.0000
 CONV       2      0      0      0     0.0500     0.0000     0.0000     0.0000
 MINI       1      0   2500      0     0.0000     0.0000     0.0000     0.0000"""

    # fixed file names
    MAE_INPUT = "input_mol.mae"
    MAE_OUTPUT = "output_mol.mae"
    SDF_OUTPUT = "output_mol.sdf"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class StepPrimeEnum:

    RECEPTOR = "receptor"  # path to the receptor MAE file
    POSEVIEWER = "poseviewer"

    # fixed file names
    SDF_INPUT = "input_mol.sdf"
    MAE_INPUT = "input_mol.mae"
    MAE_COMPLEX = "complex.mae"
    MAE_OUTPUT = "complex-out.maegz"
    MMGBSA_SCORE = "r_psp_MMGBSA_dG_Bind"
    SDF_OUTPUT = "output_mol.sdf"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class StepTurbomoleEnum:
    EXECUTION_MODE = "execution_mode"
    SUCCESS = "success"
    FAILED = "failed"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class StepPantherEnum:
    # settings.additional
    PANTHER_LOCATION = "panther_location"
    PANTHER_CONFIG_FILE = "panther_config_file"
    OUTPUT_FILE = "output_file"
    PANTHER_CONFIG_DIR = "panther_config_file"
    FIELDS = "fields"

    # fields
    FIELD_KEY_NEGATIVE_IMAGE = "negative_image"
    FIELD_KEY_COORDINATES = "5-Center"
    FIELD_KEY_PDB_FILE = "1-Pdb file"

    # parameters
    FIELDS_PARAMETERS_LIB = {
        "2-Radius": "rad.lib",
        "3-Angle": "angles.lib",
        "4-Charge": "charges.lib",
    }

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class StepKallistoEnum:

    FEATURES = "features"
    SUCCESS = "success"
    FAILURE = "failure"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class StepJazzyEnum:
    SUCCESS = "success"
    FAILURE = "failure"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class StepShaepEnum:
    # field keys for storing data
    FIELD_KEY_NEGATIVE_IMAGE = "negative_image"
    NEG_IMAGE_EXT = "mol2"

    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class StepGromacsEnum:
    FIELDS = "fields"
    COFACTORS = "cofactors"
    FORCEFIELD = "forcefield"
    INPUT_FILE = "input_file"

    FIELD_KEY_STRUCTURE = "gro"
    FIELD_KEY_PDB = "pdb"
    FIELD_KEY_TOPOL = "top"
    FIELD_KEY_TPR = "tpr"
    FIELD_KEY_MDP = "mdp"
    FIELD_KEY_XTC = "xtc"
    FIELD_KEY_ITP = "itp"
    FIELD_KEY_LOG = "log"
    FIELD_KEY_EDR = "edr"
    FIELD_KEY_NDX = "ndx"
    PROPS = "props"
    FILE_SIZE_THRESHOLD = 2000000000

    MAKE_NDX_COMMAND = "make_ndx_command"
    INDEX_FLAG = "-n"

    #  magic strings associated with ligand parametrisation step
    FORCEFIELD_ITP = "forcefield.itp"
    LIGAND_ITP = "Ligand.itp"
    INCLUDE_LIG_ITP = '#include "Ligand.itp"'
    LIG_MOLECULE_GRP = "Ligand   1\n"
    COMPLEX_TOP = "Complex.top"
    COMPLEX_PDB = "Complex.pdb"
    PROTEIN_PDB = "Protein.pdb"
    PROTEIN_TOP = "Protein.top"
    LIGAND_PDB = "Ligand.pdb"
    LIGAND_MOL2 = "Ligand.mol2"
    STD_INDEX = "index.ndx"
    STD_TOPOL = "topol.top"
    STD_TPR = "topol.tpr"
    STD_LOG = "md.log"
    STD_XTC = "traj.xtc"
    STD_TRR = "traj.trr"
    STD_STRUCTURE = "confout.gro"
    POSRE_LIG = "posre_lig.itp"
    CHARGE_METHOD = "charge_method"
    FORCE_CONSTANTS = "1000 1000 1000"
    LIG_ID = "lig_id"
    COUPLING_GROUP = "Other"
    MMPBSA_IN = "mmpbsa.in"
    THREADS = "threads"
    GROMACS_LOAD = "module load GROMACS/2021-fosscuda-2019a-PLUMED-2.7.1-Python-3.7.2"
    AMBERTOOLS_PREFIX = "ambertools_prefix"
    WATER_AND_IONS = "Water_and_ions"
    PROTEIN_OTHER = "Protein_Other"
    SIM_COMPLETE = "Finished mdrun"
    AUTO = "auto"
    TC_GRPS = "tc-grps"
    CLUSTERS_NUMBER = "clustersNumber"
    LENGTHS = "lengths"
    COUPLING_GROUPS = "coupling_groups"
    DEFAULT_MMPBSA_IN = "src/icolos/config/amber/default_mmpbsa.in"
    PARAM_METHOD = "param_method"
    MMGBSA_DG = "MMGBSA_DG"
    GAFF = "gaff"
    OPENFF = "openff"
    RESTRAINTS = "restraints"
    WATER_POSRE = """
#ifdef POSRES_WATER
[ position_restraints ]
;  i funct       fcx        fcy        fcz
   1    1       1000       1000       1000
#endif\n"""

    MULTIDIR = "multidir"
    REPLICAS = "replicas"

    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class StepOpenFFEnum:
    UNIQUE_MOLS = "unique_molecules"
    METHOD = "method"
    PARMED = "parmed"
    INTERCHANGE = "interchange"
    FORCEFIELD = "off_forcefield"


class StepCavExploreEnum:
    FIELD_KEY_DTR = "dtr"
    FIELD_KEY_CMS = "cms"
    FIELD_KEY_DX = "dx"

    # settings.additional
    CAVITY_LOCATION = "cavity_location"
    CAVITY_CONFIG_FILE = "cavity_config_file"
    OUTPUT_FILE = "output_file"
    CAVITY_CONFIG_DIR = "cavity_config_dir"
    FIELDS = "fields"
    SELECTION_TEXT = "selection_text"
    PROTEIN = "protein"
    NAME_CA = "name CA"
    FRAME_LIST_FILE = "list_of_frames.txt"
    MDPOCKET_COMMAND = "mdpocket"
    MPI_THREADS = "mpi_threads"
    EPS = "eps"
    MIN_SAMPLES = "min_samples"
    ISO_VALUE = "iso_value"
    TRAJ_TYPE = "format"
    THRESHOLD = "threshold"

    # add own fixed strings and import in step
    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class StepFepPlusEnum:
    FIELD_KEY_POSEVIEWER = "poseviewer"
    RECEPTOR_MAEGZ = "receptor.maegz"
    STRUCT_SPLIT_BASE = "split"
    STRUCTCAT_MAEGZ_OUTFILE = "concatenated.mae"
    STRUCTCAT_SDF_OUTFILE = "concatenated.sdf"
    FEP_MAPPER_OUTPUT = "out"
    FMP_OUTPUT_FILE = "out.fmp"
    LOGFILE = "multisim.log"
    EDGE_HEADER_LINE = "* Edge calculated properties (units in kcal/mol)"
    NODE_HEADER_LINE = "* Node calculated properties (units in kcal/mol)"
    SIMULATION_PROTOCOL = "* Simulation Protocol"
    SIMILARITY = "* Similarity"
    DATA_TERMINUS = "fep_mapper_cleanup: Loading output graph"
    HOST_FLAG = "-HOST"
    WAIT_FLAG = "-WAIT"
    JOBNAME_FLAG = "-JOBNAME"
    REFERENCE_DG = "ref_dg"
    JOBID_STRING = "JobId:"
    XRAY_STRUCTURES = "xray_structures"
    XRAY_SPLIT = "xray_split"
    RETRIES = "-RETRIES"

    FILE_NAME = "--name"
    FEP_EXEC_COMPLETE = "Multisim completed."
    FEP_EXEC_PARTIAL_COMPLETE = "Multisim partially completed."

    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

        # prohibit any attempt to set any values

    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class StepFilterEnum:
    FILTER_LEVEL = "filter_level"
    CONFORMERS = "conformers"
    COMPOUNDS = "compounds"
    HIGHEST_IS_BEST = "highest_is_best"
    ENUMERATIONS = "enumerations"
    CRITERIA = "criteria"
    AGGREGATION = "aggregation"
    RETURN_N = "return_n"

    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

        # prohibit any attempt to set any values

    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class StepDesmondEnum:

    PREPROCESS_MSJ = "config.msj"
    PRODUCTION_MSJ = "production.msj"
    PRODUCTION_CFG = "prod.cfg"
    MSJ_FIELDS = "msj_fields"
    CFG_FIELDS = "cfg_fields"
    SETUP_MSJ_FIELDS = "setup_msj_fields"
    CONFIG = "config"
    TOKEN_STR = "DESMOND_GPGPU:16"

    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

        # prohibit any attempt to set any values

    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class StepPdbFixerEnum:
    ADD_MISSING_HYDROGENS = "add_missing_hydrogens"
    ADD_MISSING_ATOMS = "add_missing_atoms"
    FIND_MISSING_ATOMS = "find_missing_atoms"
    FIND_MISSING_RESIDUES = "find_missing_residues"
    REPLACE_NONSTANDARD_RESIDUES = "replace_nonstandard_residues"
    REMOVE_CHAINS = "remove_chains"


class StepDSSPEnum:
    pass


class StepCressetEnum:
    SUCCESS = "success"


class StepAutoDockVinaEnum:

    ADV_RECEPTOR_PATH = "receptor_path"
    ADV_SEED = "seed"
    ADV_SEARCH_SPACE = "search_space"
    ADV_SEARCH_SPACE_CENTER_X = "--center_x"
    ADV_SEARCH_SPACE_CENTER_Y = "--center_y"
    ADV_SEARCH_SPACE_CENTER_Z = "--center_z"
    ADV_SEARCH_SPACE_SIZE_X = "--size_x"
    ADV_SEARCH_SPACE_SIZE_Y = "--size_y"
    ADV_SEARCH_SPACE_SIZE_Z = "--size_z"

    CONFIGURATION = "configuration"
    NUMBER_POSES = "number_poses"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class StepAutoDockVinaTargetPreparationEnum:

    ADV_PDBQT = ".pdbqt"
    INPUT_RECEPTOR_PDB = "input_receptor_pdb"
    OUTPUT_RECEPTOR_PDBQT = "output_receptor_pdbqt"
    PH = "pH"
    EXTRACT_BOX = "extract_box"
    EXTRACT_BOX_REFERENCE_LIGAND_PATH = "reference_ligand_path"
    EXTRACT_BOX_REFERENCE_LIGAND_FORMAT = "reference_ligand_format"
    EXTRACT_BOX_REFERENCE_LIGAND_FORMAT_PDB = "PDB"
    EXTRACT_BOX_REFERENCE_LIGAND_FORMAT_SDF = "SDF"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class StepGoldTargetPreparationEnum:

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class StepActiveLearningEnum:

    ORACLE_CONFIG = "oracle_config"
    SMILES = "SMILES"
    MOLECULE = "Moleucle"
    VIRTUAL_LIB = "virtual_lib"
    INIT_SAMPLES = "init_samples"
    MORGAN_FP = "MorganFP"
    N_ROUNDS = "n_rounds"
    BATCH_SIZE = "batch_size"
    CRITERIA = "criteria"
    VALIDATION_LIB = "validation_lib"
