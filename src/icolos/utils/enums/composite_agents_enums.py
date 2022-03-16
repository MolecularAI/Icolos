class BaseAgentEnum:

    HEADER = "header"
    STEPS = "steps"

    # header
    # ---------
    ID = "id"
    DESCRIPTION = "description"
    GLOBAL_VARIABLES = "global_variables"
    GLOBAL_SETTINGS = "global_settings"
    VERSION = "version"
    LOGGING = "logging"
    LOGGING_LOGFILE = "logfile"

    # exporting environment variables
    ENVIRONMENT = "environment"
    ENVIRONMENT_EXPORT = "export"
    ENVIRONMENT_EXPORT_KEY = "key"
    ENVIRONMENT_EXPORT_VALUE = "value"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class WorkflowEnum(BaseAgentEnum):

    WORKFLOW = "workflow"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


class SchedulerEnum(BaseAgentEnum):

    SCHEDULER = "scheduler"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")
