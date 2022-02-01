import logging

from icolos.loggers.base_logger import BaseLogger


class EntryPointLogger(BaseLogger):
    def __init__(self):
        super().__init__()

    def _initialize_logger(self):
        logger = logging.getLogger(self._LE.LOGGER_ENTRYPOINT)
        return logger
