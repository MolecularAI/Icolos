import json
import os
import numpy as np
import pandas as pd
from collections import OrderedDict
from copy import deepcopy
from typing import Tuple, List

from pydantic import BaseModel

from icolos.core.containers.compound import Conformer
from icolos.core.containers.generic import GenericData
from icolos.utils.enums.program_parameters import ModelBuilderEnum
from icolos.utils.enums.step_enums import StepModelBuilderEnum
from icolos.core.workflow_steps.io.base import StepIOBase
from icolos.core.workflow_steps.step import _LE, StepSettingsParameters
from icolos.utils.enums.write_out_enums import WriteOutEnum
from icolos.utils.execute_external.execute import Executor

_SMBE = StepModelBuilderEnum()
_SME = ModelBuilderEnum()
_WE = WriteOutEnum()


class StepModelBuilder(StepIOBase, BaseModel):
    def __init__(self, **data):
        super().__init__(**data)

        # initialize the executor
        self._initialize_backend(executor=Executor)

    def _generate_temporary_input_output_files(
        self, tmp_dir: str
    ) -> Tuple[str, str, str, str, str]:
        tmp_input_config_json = os.path.join(tmp_dir, _SMBE.TMP_INPUT_CONFIG)
        tmp_input_data_csv = os.path.join(tmp_dir, _SMBE.TMP_INPUT_DATA)
        tmp_output_best_model_pkl = os.path.join(tmp_dir, _SMBE.TMP_OUTPUT_BEST_MODEL)
        tmp_output_best_parameters_json = os.path.join(
            tmp_dir, _SMBE.TMP_OUTPUT_BEST_PARAMETERS
        )
        tmp_output_production_pkl = os.path.join(
            tmp_dir, _SMBE.TMP_OUTPUT_PRODUCTION_MODEL
        )
        return (
            tmp_input_config_json,
            tmp_input_data_csv,
            tmp_output_best_model_pkl,
            tmp_output_best_parameters_json,
            tmp_output_production_pkl,
        )

    def _update_data_block(
        self, conf: dict, tmp_input_data_csv: str, settings: StepSettingsParameters
    ) -> dict:
        # the user can specify additional things for the "data" block of the configuration
        # in the "additional" field; the input CSV file needs to be overwritten in every case, though
        specified_data_block = settings.additional.get(_SMBE.DATA, {})
        for key in specified_data_block.keys():
            conf[_SMBE.DATA][key] = specified_data_block[key]
        conf[_SMBE.DATA][_SMBE.DATA_TRAININGSET_FILE] = tmp_input_data_csv
        if _SMBE.DATA_TESTSET_FILE in conf[_SMBE.DATA].keys():
            conf[_SMBE.DATA].pop(_SMBE.DATA_TESTSET_FILE, None)
            self._logger.log(
                f"Removed test set specification, not supported yet.", _LE.WARNING
            )
        return conf

    def _write_OptunaAZ_configuration(
        self,
        tmp_input_config_json: str,
        tmp_input_data_csv: str,
        settings: StepSettingsParameters,
    ):
        config_path = settings.arguments.parameters[_SME.CONFIG]
        with open(config_path, "r") as file:
            optunaaz_conf = file.read().replace("\r", "").replace("\n", "")
            optunaaz_conf = json.loads(optunaaz_conf)
        optunaaz_conf = self._update_data_block(
            optunaaz_conf, tmp_input_data_csv, settings
        )
        with open(tmp_input_config_json, "w") as file:
            json.dump(optunaaz_conf, fp=file, indent=4)
        self._logger.log(
            f"Wrote updated OptunaAZ configuration file to {tmp_input_config_json}.",
            _LE.DEBUG,
        )

    def _write_input_csv(
        self,
        conformers: List[Conformer],
        tmp_input_data_csv: str,
        settings: StepSettingsParameters,
    ):
        def _get_tag(conformer: Conformer, tag: str) -> str:
            try:
                value = conformer.get_molecule().GetProp(tag).strip()
            except KeyError:
                value = np.nan
            return value

        smiles_column = settings.additional[_SMBE.DATA][_SMBE.DATA_INPUT_COLUMN]
        response_column = settings.additional[_SMBE.DATA][_SMBE.DATA_RESPONSE_COLUMN]

        # initialize the dictionary
        dict_result = OrderedDict()
        dict_result[_WE.RDKIT_NAME] = ["" for _ in range(len(conformers))]
        dict_result[smiles_column] = ["" for _ in range(len(conformers))]
        dict_result[response_column] = ["" for _ in range(len(conformers))]

        # populate the dictionary with the values
        for irow in range(len(conformers)):
            conf = conformers[irow]
            dict_result[_WE.RDKIT_NAME][irow] = conf.get_index_string()
            dict_result[smiles_column][irow] = _get_tag(conf, smiles_column)
            dict_result[response_column][irow] = _get_tag(conf, response_column)

        # do the writeout (after sanitation)
        df_result = pd.DataFrame.from_dict(dict_result)
        df_result.to_csv(
            path_or_buf=tmp_input_data_csv,
            sep=",",
            na_rep="",
            header=True,
            index=False,
            mode="w",
            quoting=None,
        )

    def _get_arguments(
        self,
        tmp_input_config_json: str,
        tmp_output_best_model_pkl: str,
        tmp_output_best_parameters_json: str,
        tmp_output_production_pkl: str,
    ) -> List[str]:
        arguments = [
            _SME.CONFIG,
            tmp_input_config_json,
            _SME.MERGED_MODEL_OUTPATH,
            tmp_output_production_pkl,
            _SME.BEST_MODEL_OUTPATH,
            tmp_output_best_model_pkl,
            _SME.BEST_BUILDCONFIG_OUTPATH,
            tmp_output_best_parameters_json,
        ]
        return arguments

    def _parse_output(
        self,
        tmp_input_config_json: str,
        tmp_input_data_csv: str,
        tmp_output_best_parameters_json: str,
        tmp_output_production_pkl: str,
    ):
        # loading the final model is crucial (and the end-artifact for this step)
        try:
            with open(tmp_output_production_pkl, "rb") as f:
                data = f.read()
                self.data.generic.add_file(
                    GenericData(
                        file_name=_SMBE.TMP_OUTPUT_PRODUCTION_MODEL, file_data=data
                    )
                )
        except FileNotFoundError as e:
            self._logger.log(
                f"Could not load production model from path {tmp_output_production_pkl}.",
                _LE.ERROR,
            )
            raise e

        # loading the JSON with the best hyper-parameter configuration
        try:
            with open(tmp_output_best_parameters_json, "r") as f:
                data = f.read().replace("\r", "").replace("\n", "")
                data = json.loads(data)
                self.data.generic.add_file(
                    GenericData(
                        file_name=_SMBE.TMP_OUTPUT_BEST_PARAMETERS, file_data=data
                    )
                )
        except FileNotFoundError as e:
            self._logger.log(
                f"Could not load best hyper-parameter configuration from path {tmp_output_best_parameters_json}.",
                _LE.WARNING,
            )

        # loading the input JSON for OptunaAZ
        try:
            with open(tmp_input_config_json, "r") as f:
                data = f.read()
                self.data.generic.add_file(
                    GenericData(file_name=_SMBE.TMP_INPUT_CONFIG, file_data=data)
                )
        except FileNotFoundError as e:
            self._logger.log(
                f"Could not load input CSV file from path {tmp_input_config_json}.",
                _LE.WARNING,
            )

        # loading the input CSV
        try:
            with open(tmp_input_config_json, "r") as f:
                data = f.read()
                self.data.generic.add_file(
                    GenericData(file_name=_SMBE.TMP_INPUT_DATA, file_data=data)
                )
        except FileNotFoundError as e:
            self._logger.log(
                f"Could not load input CSV file from path {tmp_input_config_json}.",
                _LE.WARNING,
            )

    def execute(self):
        # make a copy of the settings to avoid side-effects with the dictionaries
        settings = deepcopy(self.settings)

        # generate temporary files
        tmp_dir = self._move_to_temp_dir()
        (
            tmp_input_config_json,
            tmp_input_data_csv,
            tmp_output_best_model_pkl,
            tmp_output_best_parameters_json,
            tmp_output_production_pkl,
        ) = self._generate_temporary_input_output_files(tmp_dir)

        # write OptunaAZ configuration to file
        self._write_OptunaAZ_configuration(
            tmp_input_config_json=tmp_input_config_json,
            tmp_input_data_csv=tmp_input_data_csv,
            settings=settings,
        )

        # unroll all conformers
        all_conformers = []
        for compound in self.get_compounds():
            for enumeration in compound:
                all_conformers = all_conformers + enumeration.get_conformers()

        # write input CSV, derived from the conformers
        self._write_input_csv(
            conformers=all_conformers,
            tmp_input_data_csv=tmp_input_data_csv,
            settings=settings,
        )

        # execute OptunaAZ
        self._backend_executor.execute(
            command=_SME.OPTBUILD_ENTRY_POINT,
            arguments=self._get_arguments(
                tmp_input_config_json=tmp_input_config_json,
                tmp_output_best_model_pkl=tmp_output_best_model_pkl,
                tmp_output_best_parameters_json=tmp_output_best_parameters_json,
                tmp_output_production_pkl=tmp_output_production_pkl,
            ),
            check=False,
        )

        # parse the output
        self._parse_output(
            tmp_input_config_json=tmp_input_config_json,
            tmp_input_data_csv=tmp_input_data_csv,
            tmp_output_best_parameters_json=tmp_output_best_parameters_json,
            tmp_output_production_pkl=tmp_output_production_pkl,
        )

        # clean-up
        self._restore_working_dir()
        self._remove_temporary(tmp_dir)
