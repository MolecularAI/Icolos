class InputEnum:

    SOURCE_FIELD_COMPOUNDS = "compounds"
    TARGET_FIELD_COMPOUNDS = "compounds"
    TARGET_FIELD_CONFORMERS = "conformers"

    # REINVENT-compatible JSON input
    JSON_NAMES = "names"
    JSON_SMILES = "smiles"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")
