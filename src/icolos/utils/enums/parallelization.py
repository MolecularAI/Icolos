from enum import Enum


class ParallelizationEnum(str, Enum):

    STATUS_READY = "ready"
    STATUS_RUNNING = "running"
    STATUS_SUCCESS = "success"
    STATUS_FAILED = "failed"

    # try to find the internal value and return
    # def __getattr__(self, name):
    #     if name in self:
    #         return name
    #     raise AttributeError

    # # prohibit any attempt to set any values
    # def __setattr__(self, key, value):
    #     raise ValueError("No changes allowed.")
