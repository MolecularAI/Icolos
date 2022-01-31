def log_multiline_string(logger, level: str, multi_line_string: str):
    splitted = multi_line_string.split("\n")
    for line in splitted:
        logger.log(line, level)
