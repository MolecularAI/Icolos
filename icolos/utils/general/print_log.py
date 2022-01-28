import os
from icolos.loggers.blank_logger import BlankLogger


def print_log_file(path: str, logger, level):
    logger_blank = BlankLogger()
    if os.path.isfile(path):
        with open(path, "r") as log_file:
            log_file_raw = log_file.readlines()
            logger.log(f"Printing log file {path}:\n", level)
            for line in log_file_raw:
                logger_blank.log(line.rstrip("\n"), level)
            logger_blank.log("", level)
            logger.log("--- End file", level)
