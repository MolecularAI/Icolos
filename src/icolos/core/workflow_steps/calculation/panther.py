from icolos.core.containers.generic import GenericData
import os
import tempfile
import re
import numpy as np
from copy import deepcopy
from typing import List

from icolos.core.workflow_steps.calculation.base import StepCalculationBase
from icolos.utils.enums.program_parameters import PantherEnum
from icolos.utils.enums.step_enums import StepPantherEnum
from icolos.utils.execute_external.execute import Executor
from icolos.core.workflow_steps.step import _LE
from pydantic import BaseModel
from icolos.utils.general.files_paths import attach_root_path

_SPE = (
    StepPantherEnum()
)  # hold the constants to access the relevant value from initialised **data
_PE = PantherEnum()  # hold the program settings


class StepPanther(StepCalculationBase, BaseModel):

    negative_images: List = []

    def __init__(self, **data):
        super().__init__(**data)

        self._initialize_backend(executor=Executor)

    def _write_panther_config_file(self, tmp_dir):
        if not self.settings.additional[_SPE.PANTHER_CONFIG_FILE]:
            self._logger.log("No config file specified, using default.", _LE.INFO)
            panther_config = attach_root_path(
                "/icolos/config/panther/default_panther.in"
            )

        elif not os.path.isfile(self.settings.additional[_SPE.PANTHER_CONFIG_FILE]):
            self._logger.log(
                f"File not found for the provided panther config file path: {self.settings.additional[_SPE.PANTHER_CONFIG_FILE]}",
                _LE.ERROR,
            )
            raise FileNotFoundError(
                f"The specified panther config file was not found {self.settings.additional[_SPE.PANTHER_CONFIG_FILE]}"
            )

        else:
            panther_config = self.settings.additional[_SPE.PANTHER_CONFIG_FILE]

        with open(panther_config, "r") as f:
            panther_config = f.read()

        # add the parameter absolute paths to the angle, etc. file specifications
        update_dictionary = deepcopy(self.settings.additional[_SPE.FIELDS])
        update_dictionary = self._add_ligand_centroid_coordinates(update_dictionary)
        update_dictionary = self._add_parameter_locations_to_replacement_fields(
            update_dictionary
        )
        # update the configuration and write it to a file
        panther_config = self._modify_panther_config_file(
            panther_config, update_dictionary
        )

        with open(os.path.join(tmp_dir, "panther_config.in"), "w") as f:
            f.write(panther_config)

    def _add_parameter_locations_to_replacement_fields(
        self, update_dictionary: dict
    ) -> dict:
        # in case not specified (which is the main use case), use the default libraries for charges etc. that should
        # reside in the same folder as the python entry-point "panther.py"; setting absolute paths here, allows to
        # execute PANTHER for any input in any given folder
        for key, value in _SPE.FIELDS_PARAMETERS_LIB.items():
            if key not in update_dictionary.keys():
                update_dictionary[key] = os.path.join(
                    self.settings.additional[_SPE.PANTHER_LOCATION], value
                )
        return update_dictionary

    def _add_ligand_centroid_coordinates(self, update_dict: dict) -> dict:
        coordinates = self._calculate_ligand_centroid(
            self.settings.additional[_SPE.FIELDS][_SPE.FIELD_KEY_PDB_FILE]
        )
        update_dict[_SPE.FIELD_KEY_COORDINATES] = coordinates
        return update_dict

    def _calculate_ligand_centroid(self, file: str) -> str:
        with open(file, "r") as f:
            file_lines = f.readlines()
        file_lines = [
            line for line in file_lines if "X   0" in line and len(line.split()) > 5
        ]
        if not file_lines:
            self._logger.log(
                "No lines corresponding to the ligand found! Centroid will not be correct",
                _LE.WARNING,
            )
        a = np.genfromtxt(file_lines, usecols=[6, 7, 8], skip_header=1)
        avg = list(a.mean(axis=0))
        avg = [str(i) for i in avg]
        return " ".join(avg)

    def _modify_panther_config_file(
        self, config_file: str, update_dictionary: dict
    ) -> str:
        for key, value in update_dictionary.items():
            pattern = rf"({key}.*:: ).*"
            pattern = re.compile(pattern)
            config_file = re.sub(pattern, rf"\1 {value}", config_file)
        return config_file

    def _execute_backend(self, tmp_dir):
        arguments = [
            os.path.join(
                self.settings.additional[_SPE.PANTHER_LOCATION], _PE.PANTHER_ENTRYPOINT
            ),
            os.path.join(tmp_dir, _PE.PANTHER_CONFIG),
            os.path.join(tmp_dir, _PE.PANTHER_OUTPUT_FILE),
        ]
        self._backend_executor.execute(
            command=_PE.PANTHER_PTYHON2, arguments=arguments, check=True
        )

    def _parse_panther_output(self, tmp_dir):
        try:
            with open(os.path.join(tmp_dir, _PE.PANTHER_OUTPUT_FILE), "r") as f:
                data = f.read()
                self.data.generic.add_file(
                    GenericData(file_name=_PE.PANTHER_OUTPUT_FILE, file_data=data)
                )
        except FileNotFoundError:
            self._logger.log(
                f"No panther output file was produced for step {self.step_id}, subsequent steps that depend on the negative image will fail.",
                _LE.WARNING,
            )

    def execute(self):
        tmp_dir = self._make_tmpdir()
        self._write_panther_config_file(tmp_dir)
        self._execute_backend(tmp_dir)
        self._logger.log("Executed PANTHER and obtained negative image.", _LE.INFO)
        self._logger.log(
            f"Calculated negative image for configuration file in {tmp_dir}.", _LE.DEBUG
        )
        self._parse_panther_output(tmp_dir)
        self._remove_temporary(tmp_dir)
