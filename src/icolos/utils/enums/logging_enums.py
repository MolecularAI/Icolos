class LoggingConfigEnum:

    # set levels (for now, they match to the "logging" default ones)
    DEBUG = "debug"
    INFO = "info"
    WARNING = "warning"
    ERROR = "error"
    EXCEPTION = "exception"

    # paths to the configuration JSONs that are shipped with Icolos
    PATH_CONFIG_DEFAULT = "src/icolos/config/logging/default.json"
    PATH_CONFIG_VERBOSE = "src/icolos/config/logging/verbose.json"
    PATH_CONFIG_DEBUG = "src/icolos/config/logging/debug.json"
    PATH_CONFIG_TUTORIAL = "src/icolos/config/logging/tutorial.json"

    # high-level loggers defined in the configurations
    LOGGER_IO = "io"
    LOGGER_STEP = "step"
    LOGGER_AGENT = "agent"
    LOGGER_ENTRYPOINT = "entrypoint"
    LOGGER_BLANK = "blank"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")
