def get_progress_bar_string(
    done, total, prefix="", suffix="", decimals=1, length=100, fill="█"
):
    percent = ("{0:." + str(decimals) + "f}").format(100 * (done / float(total)))
    filledLength = int(length * done // total)
    bar = fill * filledLength + "-" * (length - filledLength)
    return f"{prefix}|{bar}| {percent}% {suffix}"
