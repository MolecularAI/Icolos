class WriteOutEnum:

    RDKIT_NAME = "_Name"
    INDEX_STRING = "index_string"
    COMPOUND_NAME = "compound_name"

    # REINVENT-compatible JSON write-out
    JSON_RESULTS = "results"
    JSON_NAMES = "names"
    JSON_NA = ""
    JSON_VALUES = "values"
    JSON_VALUES_KEY = "values_key"

    SDF = "sdf"
    PDB = "pdb"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class RunVariablesEnum:

    PREFIX = "["
    POSTFIX = "]"
    COMPOUND_ID = "compound_id"
    ENUMERATION_ID = "enumeration_id"
    CONFORMER_ID = "conformer_id"
    COMPOUND_NAME = "compound_name"
    ENUMERATION_STRING = "enumeration_string"
    CONFORMER_STRING = "conformer_string"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")
