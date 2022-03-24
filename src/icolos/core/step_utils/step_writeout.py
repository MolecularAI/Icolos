import os
from collections import OrderedDict
from copy import deepcopy

import numpy as np
import pandas as pd
import json
from typing import List
from pydantic import BaseModel, PrivateAttr
from rdkit import Chem
from pathlib import Path

# from icolos.core.composite_agents.workflow import WorkflowData

from icolos.core.containers.compound import Compound, Conformer
from icolos.core.step_utils.input_preparator import StepData
from icolos.core.step_utils.run_variables_resolver import RunVariablesResolver
from icolos.loggers.steplogger import StepLogger
from icolos.utils.enums.logging_enums import LoggingConfigEnum
from icolos.utils.enums.step_enums import StepBaseEnum, StepGromacsEnum
from icolos.utils.enums.write_out_enums import WriteOutEnum

_WE = WriteOutEnum()
_LE = LoggingConfigEnum()
_SBE = StepBaseEnum
_SGE = StepGromacsEnum()


class StepWriteoutCompoundAggregationParameters(BaseModel):
    mode: _SBE = _SBE.WRITEOUT_COMP_AGGREGATION_MODE_ALL
    highest_is_best: bool = True
    key: str = None


class StepWriteoutCompoundParameters(BaseModel):
    category: _SBE
    aggregation: StepWriteoutCompoundAggregationParameters = (
        StepWriteoutCompoundAggregationParameters()
    )
    key: str = None
    selected_tags: List[str] = None


class StepWriteoutGenericParameters(BaseModel):
    key: str


class StepWriteoutGromacsParameters(BaseModel):
    key: str


class StepWriteoutDestinationParameters(BaseModel):
    resource: str = None
    type: _SBE = _SBE.WRITEOUT_DESTINATION_TYPE_FILE
    format: _SBE = _SBE.FORMAT_TXT
    merge: bool = True
    mode: _SBE = _SBE.WRITEOUT_DESTINATION_BASE_NAME


class StepWriteoutParameters(BaseModel):
    compounds: StepWriteoutCompoundParameters = None
    generic: StepWriteoutGenericParameters = None
    gmx_state: StepWriteoutGromacsParameters = None
    destination: StepWriteoutDestinationParameters = None


class WriteOutHandler(BaseModel):

    config: StepWriteoutParameters
    data: StepData = None
    workflow_data: BaseModel = None

    class Config:
        underscore_attrs_are_private = True

    _logger = PrivateAttr()

    def __init__(self, **data):
        super().__init__(**data)
        self._logger = StepLogger()

    def set_data(self, data: StepData):
        self.data = deepcopy(data)

    def set_workflow_data(self, data):
        self.workflow_data = data

    def get_data(self) -> StepData:
        return self.data

    def _handle_destination_type(self):
        if self.config.destination.type.lower() in (
            _SBE.WRITEOUT_DESTINATION_TYPE_FILE,
            _SBE.WRITEOUT_DESTINATION_TYPE_REINVENT,
            _SBE.WRITEOUT_DESTINATION_DIR,
        ):
            return self.config.destination.resource
        elif (
            self.config.destination.type.lower() == _SBE.WRITEOUT_DESTINATION_TYPE_REST
        ):
            raise ValueError("REST end-point destination type not supported yet.")
        raise ValueError(
            f"Destination type {self.config.destination.type} not supported."
        )

    def _write_compounds(self):
        resource = self._handle_destination_type()
        resolver = RunVariablesResolver()
        if self.config.compounds.category == _SBE.WRITEOUT_COMP_CATEGORY_CONFORMERS:
            if self.config.destination.format.upper() == _SBE.FORMAT_CSV:
                if self.config.destination.format.upper() != _SBE.FORMAT_CSV:
                    raise NotImplementedError(
                        "Only supporting CSV write-out format for tabular data."
                    )
                self._writeout_tabular()
            elif self.config.destination.format.upper() == _SBE.FORMAT_JSON:
                self._writeout_reinvent()
            elif self.config.destination.format.upper() == _SBE.FORMAT_SDF:

                def _write_compounds(compounds: List[Compound], resource: str):
                    # TODO: deal with resolving resources differently (also for writing enumerations below)
                    resource_resolved = resource
                    for compound in compounds:
                        for enum in compound.get_enumerations():
                            if len(enum.get_conformers()) > 0:
                                resource_resolved = resolver.resolve(resource, enum[0])
                                break
                    self._make_folder(resource_resolved)

                    writer = Chem.SDWriter(resource_resolved)
                    written = 0

                    for comp in compounds:
                        for enum in comp:
                            for conf in enum:
                                molecule = conf.get_molecule()
                                if (
                                    comp.get_name() is not None
                                    and comp.get_name() != ""
                                ):
                                    molecule.SetProp(_WE.COMPOUND_NAME, comp.get_name())
                                molecule.SetProp(
                                    _WE.RDKIT_NAME, conf.get_index_string()
                                )
                                molecule.SetProp(
                                    _WE.INDEX_STRING, conf.get_index_string()
                                )
                                writer.write(molecule)
                                written += 1
                    writer.close()
                    self._logger.log(
                        f"Wrote {written} conformers to file {resource_resolved}.",
                        _LE.DEBUG,
                    )

                # TODO: At the moment, this only splits at the compound level (taking the first conformer for resolving),
                if self.config.destination.merge:
                    _write_compounds(self.data.compounds, resource=resource)
                else:
                    for comp in self.data.compounds:
                        _write_compounds([comp], resource)
        elif self.config.compounds.category == _SBE.WRITEOUT_COMP_CATEGORY_ENUMERATIONS:
            if not self.config.destination.format.upper() == _SBE.FORMAT_SDF:
                raise NotImplementedError(
                    "This write-out is not supported for enumerations."
                )
            else:

                def _write_compounds(compounds: List[Compound], resource: str):
                    # TODO: deal with resolving resources differently (also for writing conformers above)
                    resource_resolved = resource
                    for compound in compounds:
                        if len(compound.get_enumerations()) > 0:
                            resource_resolved = resolver.resolve(resource, compounds[0])
                            break

                    self._make_folder(resource_resolved)
                    writer = Chem.SDWriter(resource_resolved)
                    written = 0
                    for comp in compounds:
                        for enum in comp:
                            molecule = enum.get_molecule()
                            if comp.get_name() is not None and comp.get_name() != "":
                                molecule.SetProp(_WE.COMPOUND_NAME, comp.get_name())
                            molecule.SetProp(_WE.RDKIT_NAME, enum.get_index_string())
                            molecule.SetProp(_WE.INDEX_STRING, enum.get_index_string())
                            writer.write(molecule)
                            written += 1
                    writer.close()
                    self._logger.log(
                        f"Wrote {written} enumeration molecules to file {resource_resolved}.",
                        _LE.DEBUG,
                    )

                if self.config.destination.merge:
                    _write_compounds(self.data.compounds, resource=resource)
                else:
                    for comp in self.data.compounds:
                        _write_compounds([comp], resource)
        elif self.config.compounds.category == _SBE.WRITEOUT_COMP_CATEGORY_EXTRADATA:
            if self.config.destination.format.upper() != _SBE.FORMAT_TXT:
                raise ValueError(
                    f"For writing out extra-data (attached to conformers), only TXT is supported as format."
                )
            # TODO: Does merging here makes any sense?
            for comp in self.data.compounds:
                for enum in comp:
                    for conf in enum:
                        resource_resolved = resolver.resolve(resource, conf)
                        self._make_folder(resource_resolved)
                        with open(resource_resolved, "w") as f:
                            content = conf.get_extra_data()[self.config.compounds.key]
                            if isinstance(content, list):
                                for line in content:
                                    f.write(line.rstrip("\n") + "\n")
                            elif isinstance(content, str):
                                f.write(content)
                            else:
                                raise ValueError(
                                    "Extra data must be either a string or a list of strings."
                                )
        else:
            raise ValueError(f"{self.config.compounds.category} not supported.")

    def _write_generic_data(self):
        # type and format do not apply here, simply overwrite defaults
        self.config.destination.type = _SBE.WRITEOUT_DESTINATION_TYPE_FILE
        self.config.destination.format = _SBE.FORMAT_TXT
        resource = self._handle_destination_type()
        self._make_folder(resource)
        if self.config.destination.mode == _SBE.WRITEOUT_DESTINATION_DIR:
            # The output path should be a directory only
            assert not os.path.isfile(resource)
            os.makedirs(resource, exist_ok=True)
        # write out all files from that step with the required extension
        for idx, file in enumerate(
            self.data.generic.get_files_by_extension(self.config.generic.key)
        ):
            if self.config.destination.mode == _SBE.WRITEOUT_DESTINATION_BASE_NAME:
                parts = resource.split(".")
                resource = parts[0] + f"_{idx}." + parts[1]
                file.write(resource, join=False)
            elif self.config.destination.mode == _SBE.WRITEOUT_DESTINATION_AUTOMATIC:
                # take the original file name from the step (these tend not to be very descriptive)
                parts = file.get_file_name().split(".")
                file_name = parts[0] + f"_{idx}." + parts[1]
                resource = os.path.join("/".join(resource.split("/")[:-1]), file_name)
                file.write(resource, join=False)
            elif self.config.destination.mode == _SBE.WRITEOUT_DESTINATION_DIR:
                resource = resource
                assert os.path.isdir(resource)
                file.write(resource, join=True, final_writeout=True)

    def _write_gromacs_data(self):
        """
        Handle writeout from gromacs topology state
        """
        self.config.destination.type = _SBE.WRITEOUT_DESTINATION_TYPE_FILE
        self.config.destination.format = _SBE.FORMAT_TXT
        self.config.destination.type = _SBE.WRITEOUT_DESTINATION_DIR
        resource = self._handle_destination_type()
        os.makedirs(resource, exist_ok=True)
        writeout_keys = map(lambda s: s.strip(), self.config.gmx_state.key.split(","))
        for key in writeout_keys:
            if key == _SGE.FIELD_KEY_TOPOL:
                self.data.gmx_state.write_topol(resource)
            elif key == _SGE.FIELD_KEY_NDX:
                self.data.gmx_state.write_ndx(resource)
            elif key == _SGE.PROPS:
                self.data.gmx_state.write_props(resource)
            elif key == _SGE.FIELD_KEY_LOG:
                self.data.gmx_state.write_log(resource)
            elif key == _SGE.FIELD_KEY_TPR:
                self.data.gmx_state.write_tpr(resource)
            elif key == _SGE.FIELD_KEY_XTC:
                # if we have multiple trajectories, write them out sequentially, with index attached
                if len(self.data.gmx_state.trajectories.keys()) > 1:
                    for k, v in self.data.gmx_state.trajectories.items():
                        parts = v.get_file_name().split(".")
                        file_name = parts[0] + "_" + str(k) + "." + parts[1]
                        self.data.gmx_state.write_trajectory(
                            resource, file=file_name, index=k
                        )

                else:
                    self.data.gmx_state.write_trajectory(resource)
            elif key == _SGE.FIELD_KEY_STRUCTURE:
                if len(self.data.gmx_state.structures.keys()) > 1:
                    for k, v in self.data.gmx_state.structures.items():
                        parts = v.get_file_name().split(".")
                        file_name = parts[0] + "_" + str(k) + "." + parts[1]
                        self.data.gmx_state.write_structure(
                            resource, file=file_name, index=k
                        )

                else:
                    self.data.gmx_state.write_structure(resource)
            else:
                raise ValueError(
                    f"Gromacs file of type {key} is not supported for writeout"
                )

    def write(self):
        if (
            self.config.compounds is not None
            and self.config.generic is not None
            and self.config.gromacs_state is not None
        ):
            raise ValueError("Only specify one type of writeout per block!")

        if self.config.compounds is not None:
            self._write_compounds()
        elif self.config.generic is not None:
            self._write_generic_data()
        elif self.config.gmx_state is not None:
            self._write_gromacs_data()
        else:
            raise ValueError(
                "Either compounds, generic or gromacs data has to be specified."
            )

    def _writeout_reinvent(self):
        def _get_conf_by_comp_name(confs: List[Conformer], comp_name: str) -> Conformer:
            # assumes there is at most 1 conformer / compound left at this stage, as is required by REINVENT
            for conf in confs:
                if conf.get_compound_name() == comp_name:
                    return conf
            return None

        dict_result = {_WE.JSON_RESULTS: []}
        tags = self._get_selected_tags()

        # add names, including those for which no conformer has been obtained
        dict_result[_WE.JSON_NAMES] = [comp.get_name() for comp in self.data.compounds]

        # do aggregation (might remove conformers)
        confs_unrolled = self._apply_aggregation(self.data.compounds)

        # add values (derived from molecule tags)
        # TODO: if no conformers are left, we need to write out an empty JSON that tells REINVENT that none worked
        for tag in tags:
            values = []
            for comp_name in dict_result[_WE.JSON_NAMES]:
                conf = _get_conf_by_comp_name(confs=confs_unrolled, comp_name=comp_name)
                if conf is not None:
                    try:
                        value = conf.get_molecule().GetProp(tag)
                    except KeyError:
                        value = _WE.JSON_NA
                else:
                    value = _WE.JSON_NA
                values.append(value.strip())
            dict_result[_WE.JSON_RESULTS].append(
                {_WE.JSON_VALUES_KEY: tag, _WE.JSON_VALUES: values}
            )

        # TODO: refactor that part
        resource = self._handle_destination_type()
        if len(confs_unrolled) > 0:
            resolver = RunVariablesResolver()
            resource_resolved = resolver.resolve(resource, confs_unrolled[0])
        else:
            resource_resolved = resource
            self._logger.log(
                f"No conformers obtained, write-out resource resolving disabled.",
                _LE.WARNING,
            )
        self._make_folder(resource_resolved)

        # write-out according to destination type
        # TODO: there seems to be an issue here, when multiple write-out blocks are specified and no conformers are
        #       left: only the first block gets executed and if that's not the REINVENT one, the run will crash
        if self.config.destination.type.lower() in (
            _SBE.WRITEOUT_DESTINATION_TYPE_REINVENT,
            _SBE.WRITEOUT_DESTINATION_TYPE_FILE,
        ):
            with open(resource_resolved, "w") as f:
                json.dump(dict_result, f, indent=4)
        elif self.config.destination.type.lower() in (
            _SBE.WRITEOUT_DESTINATION_TYPE_STDOUT,
            _SBE.WRITEOUT_DESTINATION_TYPE_STDERR,
        ):
            json.dump(dict_result, resource_resolved, indent=4)
        else:
            raise ValueError(
                f"Destination type {self.config.destination.type} not supported for this function."
            )

    def _get_selected_tags(self) -> List[str]:
        # this function returns a list of tags (strings) that are to be considered for e.g. tabular write-out
        # if the respective configuration field is set to "None", use all tags (over all compounds in a batch)
        if self.config.compounds.selected_tags is not None:
            if isinstance(self.config.compounds.selected_tags, list):
                list_tags = self.config.compounds.selected_tags
            elif isinstance(self.config.compounds.selected_tags, str):
                list_tags = [self.config.compounds.selected_tags]
            else:
                raise ValueError(
                    f'Tag selection "{self.config.compounds.selected_tags}" set to illegal value.'
                )
        else:
            # get all tags for all compounds
            list_tags = []
            for comp in self.data.compounds:
                for enum in comp:
                    for conf in enum:
                        list_tags = list_tags + list(conf.get_molecule().GetPropNames())

        list_tags = list(set(list_tags))
        return list_tags

    def _initialize_dict_csv(
        self, keys: List[str], nrow: int, fill_value=np.NaN
    ) -> OrderedDict:
        return_dict = OrderedDict()
        for key in keys:
            return_dict[key] = [fill_value for _ in range(nrow)]
        return return_dict

    def _apply_aggregation(self, compounds: List[Compound]) -> List[Conformer]:
        if (
            self.config.compounds.aggregation.mode
            == _SBE.WRITEOUT_COMP_AGGREGATION_MODE_ALL
        ):
            return self._unroll_conformers(compounds)

        confs_remaining = []
        if (
            self.config.compounds.aggregation.mode
            == _SBE.WRITEOUT_COMP_AGGREGATION_MODE_BESTPERENUMERATION
        ):
            raise NotImplementedError("Best per enumeration is not yet implemented.")
        elif (
            self.config.compounds.aggregation.mode
            == _SBE.WRITEOUT_COMP_AGGREGATION_MODE_BESTPERCOMPOUND
        ):
            for comp in compounds:
                unrolled_conformers = self._unroll_conformers([comp])
                if len(unrolled_conformers) == 0:
                    continue
                values = [
                    float(
                        conf.get_molecule().GetProp(
                            self.config.compounds.aggregation.key
                        )
                    )
                    for conf in unrolled_conformers
                ]
                index_best = (
                    values.index(max(values))
                    if self.config.compounds.aggregation.highest_is_best
                    else values.index(min(values))
                )
                confs_remaining.append(unrolled_conformers[index_best])
        return confs_remaining

    def _unroll_conformers(self, compounds: List[Compound]) -> List[Conformer]:
        result = []
        for comp in compounds:
            for enum in comp:
                for conf in enum:
                    result.append(conf)
        return result

    def _writeout_tabular(self):
        # get all tags of the molecules that are to be considered
        tags = self._get_selected_tags()

        # remove the compound_name and _Name, as they will be specifically added at the beginning
        if _WE.COMPOUND_NAME in tags:
            tags.remove(_WE.COMPOUND_NAME)
        if _WE.RDKIT_NAME in tags:
            tags.remove(_WE.RDKIT_NAME)

        # do aggregation (might remove conformers)
        confs_unrolled = self._apply_aggregation(self.data.compounds)

        # initialize a dictionary with all tags as keys and filled with NA for every position
        dict_result = self._initialize_dict_csv(
            keys=[_WE.RDKIT_NAME, _WE.COMPOUND_NAME] + tags, nrow=len(confs_unrolled)
        )

        # resolve resource
        # TODO: refactor that part
        resource = self._handle_destination_type()
        resolver = RunVariablesResolver()
        if len(confs_unrolled) == 0:
            raise ValueError("No conformers found.")
        resource_resolved = resolver.resolve(resource, confs_unrolled[0])
        self._make_folder(resource_resolved)

        # populate the dictionary with the values (if present)
        for irow in range(len(confs_unrolled)):
            # add the internal Icolos identifier
            conf = confs_unrolled[irow]
            dict_result[_WE.RDKIT_NAME][irow] = conf.get_index_string()

            # add the compound name, if specified
            name = conf.get_compound_name()
            dict_result[_WE.COMPOUND_NAME][irow] = "" if name is None else name
            for tag in tags:
                try:
                    value = conf.get_molecule().GetProp(tag).strip()
                except KeyError:
                    value = np.nan
                dict_result[tag][irow] = value

        # do the writeout (after sanitation)
        df_result = pd.DataFrame.from_dict(dict_result)
        df_result = self._sanitize_df_columns(df=df_result)
        df_result.to_csv(
            path_or_buf=resource_resolved,
            sep=",",
            na_rep="",
            header=True,
            index=False,
            mode="w",
            quoting=None,
        )
        self._logger.log(
            f"Wrote data frame with {len(confs_unrolled)} rows and {len(tags)} columns to file {resource_resolved}.",
            _LE.DEBUG,
        )

    def _sanitize_df_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        cols_before = df.columns.to_list()
        df.columns = (
            df.columns.str.strip()
            .str.replace(" ", "_")
            .str.replace("(", "")
            .str.replace(")", "")
            .str.replace("/", "_")
            .str.replace("[", "")
            .str.replace("]", "")
        )
        for col_before, col_after in zip(cols_before, df.columns.to_list()):
            if col_before != col_after:
                self._logger.log(
                    f"Sanitized column name {col_before} to {col_after}.", _LE.WARNING
                )
        return df

    def _make_folder(self, path):
        if isinstance(path, str):
            if not os.path.isdir(path):
                path = os.path.dirname(path)
            Path(path).mkdir(parents=True, exist_ok=True)
