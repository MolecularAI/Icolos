class ExecutorEnum:

    RUNTIME_GLOBAL_VARIABLE_WORKDIR = "work_dir"
    RUNTIME_GLOBAL_VARIABLE_ENTRYPOINTDIR = "entrypoint_dir"
    RUNTIME_GLOBAL_VARIABLE_CONFIGDIR = "config_dir"
    RUNTIME_GLOBAL_VARIABLE_PACKAGEDIR = "package_dir"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")
