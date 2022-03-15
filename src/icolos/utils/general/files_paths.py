import os
import shutil
import time
import json
import tempfile
from typing import Tuple

from icolos.utils.enums.step_enums import StepBaseEnum
from icolos.utils.enums.general_utils_enums import CheckFileGenerationEnum

_SE = StepBaseEnum
_FG = CheckFileGenerationEnum()


def check_file_availability(
    path: str, interval_sec: int = 1, maximum_sec: int = 10
) -> str:
    counter = 0
    while not os.path.exists(path):
        # wait for an interval
        time.sleep(interval_sec)
        counter = counter + 1

        # if there's time left, proceed
        if maximum_sec is not None and (counter * interval_sec) > maximum_sec:
            break
    if os.path.exists(path):
        if os.path.getsize(path) == 0:
            return _FG.GENERATED_EMPTY
        else:
            return _FG.GENERATED_SUCCESS
    else:
        return _FG.NOT_GENERATED


def remove_folder(folder_path: str):
    if os.path.isdir(folder_path):
        shutil.rmtree(folder_path)


def empty_output_dir(path: str):
    for root, subf, files in os.walk(path):
        for file in files:
            os.remove(os.path.join(root, file))


def move_up_directory(path, n=1):
    """Function, to move up 'n' directories for a given "path"."""
    # add +1 to take file into account
    if os.path.isfile(path):
        n += 1
    for _ in range(n):
        path = os.path.dirname(os.path.abspath(path))
    return path


def attach_root_path(path, n=4):
    """Function to attach the root path of the module for a given "path"."""
    ROOT_DIR = move_up_directory(os.path.abspath(__file__), n=n)
    return os.path.join(ROOT_DIR, path)


def lines_in_file(path):
    with open(path) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def dict_from_json_file(path):
    with open(path, "r") as f:
        return json.load(f)


def any_in_file(path, strings):
    if isinstance(strings, str):
        strings = [strings]
    if os.path.isfile(path):
        with open(path, "r") as f:
            file_raw = f.readlines()
            for string in strings:
                if any(string in line for line in file_raw):
                    return True
            return False
    else:
        return False


def infer_input_type(path: str) -> str:
    basename = os.path.basename(path)
    ending = basename[-3:].upper()
    if ending in [_SE.FORMAT_SDF, _SE.FORMAT_CSV, _SE.FORMAT_SMI]:
        return ending
    else:
        raise ValueError(f"Ending {ending} not supported.")


def gen_tmp_file(
    suffix: str = None, prefix: str = None, dir: str = None, text: bool = True
) -> Tuple[str, str]:
    """Function wraps tempfile.mkstemp(), but closes the connection and returns the file name instead of the handler."""
    # note that in contrast to the underlying "mkstemp" function, "text" is set to True here
    fhandle, path = tempfile.mkstemp(suffix=suffix, prefix=prefix, dir=dir, text=text)
    os.close(fhandle)
    return os.path.basename(path), path
