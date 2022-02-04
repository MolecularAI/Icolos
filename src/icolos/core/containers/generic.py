from shutil import copyfile
from distutils.dir_util import copy_tree
import json
import os
import sys
from typing import Any, List, Dict, Union
from copy import Error


class GenericData:
    """Container class to hold generic data of any file type"""

    def __init__(
        self,
        file_name: str,
        file_data=None,
        argument=True,
        file_id: int = None,
        extension: str = None,
    ):
        self._extension = (
            extension if extension is not None else file_name.split(".")[-1]
        )
        self._file_name = file_name
        self._file_data = file_data
        self._file_id = file_id
        # self._argument: bool = argument
        self._file_size = self.calculate_file_size()

    def get_file_name(self) -> str:
        return self._file_name

    def get_data(self) -> Any:
        return self._file_data

    def calculate_file_size(self):
        return sys.getsizeof(self._file_data)

    def get_extension(self):
        return self._extension

    def set_data(self, data):
        self._file_data = data

    def set_file_name(self, file_name):
        self._file_name = file_name

    def set_id(self, file_id):
        self._file_id = file_id

    def get_id(self):
        return self._file_id

    def set_extension(self, extension):
        self._extension = extension

    def write(self, path: str, join: bool = True, final_writeout: bool = False):
        """
        Handles all I/O operations for generic data.  Support for handling directories and symlinks
        """
        orig_path = path
        if join:
            path = os.path.join(path, self.get_file_name())

        if str(self._file_data).startswith("/"):
            # file data is a path, copy the file to the destination
            # if it's a file, its stored like this because it's large (> 2GB)
            if os.path.isfile(self._file_data):
                if not final_writeout:
                    # if this is a writeout to a step, we can simply create a simlink
                    os.symlink(self._file_data, path, target_is_directory=False)
                else:
                    # we cannot do this for the final writeout since /scratch or /tmp will eventually get cleaned
                    copyfile(self._file_data, path)

            elif os.path.isdir(self._file_data):
                # copy the entire directory to the parent dir
                copy_tree(self._file_data, orig_path)
        elif isinstance(self._file_data, list):
            with open(path, "w") as f:
                f.writelines(self._file_data)

        elif isinstance(self._file_data, str):
            with open(path, "w") as f:
                f.write(self._file_data)
        elif isinstance(self._file_data, dict):
            with open(path, "w") as f:
                f.write(json.dumps(self._file_data))
        else:
            with open(path, "wb") as f:
                f.write(self._file_data)

    def update_data(self, data):
        if sys.getsizeof(data) != self._file_size:
            self.set_data(data)

    def __repr__(self):
        return f"GenericData object - name: {self._file_name}, extension: {self._extension}."

    def __str__(self):
        return self.__repr__()


class GenericContainer:
    """Container class to hold the instances of the Generic class, separated by extension"""

    def __init__(self):
        self._file_dict: Dict[str, List] = {}

    # self._paths = []
    # self._strings = []

    def add_file(self, file: GenericData):
        ext = file.get_extension()
        file.set_id(self.get_next_file_id(ext))
        try:
            self._file_dict[ext].append(file)
        except NameError:
            self._initialise_list(ext)
            self._file_dict[ext].append(file)

    def _initialise_list(self, ext):
        self._file_dict[ext] = []

    def get_next_file_id(self, ext):
        ids = [file.get_id() for file in self.get_files_by_extension(ext)]
        if len(ids) == 0:
            return 0
        else:
            return max(ids) + 1

    def get_file_by_index(self, index):
        for file in self.get_flattened_files():
            if file.get_id() == index:
                return file

    def add_files(self, files: List[GenericData]):
        extensions = list(set([f.get_extension() for f in files]))
        if len(extensions) > 1:
            raise Error("Cannot have more than one type of file")
        else:
            if extensions[0] in self._file_dict.keys():
                self._file_dict[extensions[0]].extend(files)
            else:
                self._file_dict[extensions[0]] = files

    def get_all_files(self) -> Dict[str, List]:
        return self._file_dict

    def get_files_by_extension(self, ext: str) -> List[GenericData]:
        if ext in self._file_dict.keys():
            return self._file_dict[ext]
        else:
            self._initialise_list(ext)
            return self._file_dict[ext]

    def get_file_names_by_extension(self, ext: str):
        try:
            return [f.get_file_name() for f in self._file_dict[ext]]
        except KeyError:
            self._initialise_list(ext)
            return [f.get_file_name() for f in self._file_dict[ext]]

    def get_file_types(self):
        return self._file_dict.keys()

    def get_flattened_files(self) -> List[GenericData]:
        rtn_files = []
        for key in self._file_dict.keys():
            for file in self._file_dict[key]:
                rtn_files.append(file)
        return rtn_files

    def get_file_by_name(self, name):
        for file in self.get_flattened_files():
            if file.get_file_name() == name:
                return file

    def clear_file_dict(self):
        self._file_dict = {}

    def get_argument_by_extension(
        self, ext, rtn_file_object=False
    ) -> Union[GenericData, str]:
        files = []
        for file in self.get_flattened_files():
            if file.get_extension() == ext:
                files.append(file)
        assert len(files) > 0, f"No files with extension {ext} were found!"
        try:
            assert len(files) == 1
        except AssertionError:
            print(
                f"Found multiple files with extension {ext}, select the index of the file to be passed as an argument\n"
            )
            print("######################")
            for idx, file in enumerate(files):
                print(f"{idx}: {file.get_file_name()}")
            print("######################")
            index = input(">>> ")
            files = [files[int(index)]]

        if not rtn_file_object:
            return files[0].get_file_name()
        else:
            return files[0]

    def write_out_all_files(self, folder):
        """flattens all files in the container and writes to the specified directory"""
        for file in self.get_flattened_files():
            file.write(folder)
