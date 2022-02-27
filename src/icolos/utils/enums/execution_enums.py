from enum import Enum


class ExecutionPlatformEnum(str, Enum):
    LOCAL = "local"
    SLURM = "slurm"
    PARTITION = "partition"
    TIME = "time"
    GRES = "gres"
    MEM = "mem"
    CORES = "cores"
    CORE = "core"
    GPU = "gpu"
    SERIAL = "serial"
