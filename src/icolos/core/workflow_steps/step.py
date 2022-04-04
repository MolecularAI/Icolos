from subprocess import CompletedProcess
import time

from icolos.core.containers.generic import GenericContainer, GenericData
import multiprocessing
import shutil
import tempfile
from typing import Callable, List, Dict, Tuple

from pydantic import BaseModel, PrivateAttr
from rdkit import Chem
from copy import deepcopy
import os
from icolos.core.containers.gmx_state import GromacsState
from icolos.core.containers.perturbation_map import PerturbationMap


from icolos.core.step_utils.input_preparator import (
    StepData,
    InputPreparator,
    StepInputParameters,
)
from icolos.loggers.steplogger import StepLogger
from icolos.loggers.blank_logger import BlankLogger
from icolos.utils.enums.step_enums import StepGromacsEnum
from icolos.core.containers.compound import Compound, Conformer
from icolos.core.step_utils.step_writeout import (
    StepWriteoutParameters,
    WriteOutHandler,
    _SBE,
)
from icolos.utils.enums.execution_enums import (
    ExecutionPlatformEnum,
)
from icolos.utils.execute_external.execute import Executor
from icolos.utils.execute_external.slurm_executor import SlurmExecutor
from icolos.utils.general.icolos_exceptions import StepFailed

from icolos.utils.enums.compound_enums import CompoundTagsEnum
from icolos.utils.enums.logging_enums import LoggingConfigEnum
from icolos.utils.enums.write_out_enums import WriteOutEnum
from icolos.utils.general.files_paths import gen_tmp_file, any_in_file
from icolos.utils.general.parallelization import SubtaskContainer, Subtask
from tempfile import mkdtemp
from distutils.dir_util import copy_tree
from icolos.core.containers.compound import unroll_enumerations, unroll_conformers
from icolos.utils.general.progress_bar import get_progress_bar_string

_LE = LoggingConfigEnum()
_WE = WriteOutEnum()
_CTE = CompoundTagsEnum()
_SGE = StepGromacsEnum()
_EPE = ExecutionPlatformEnum


class StepFailurePolicyParameters(BaseModel):
    n_tries: int = 1
    retry_wait_seconds: int = 10


class StepExecutionResourceParameters(BaseModel):
    partition: str = _EPE.CORE
    time: str = "12:00:00"
    gres: str = None
    tasks: str = None
    mem: str = None
    cores: int = None
    modules: List = []
    other_args: dict = {}
    additional_lines: List = []


class StepExecutionParameters(BaseModel):
    class StepExecutionParallelizationParameters(BaseModel):
        jobs: int = 1
        max_length_sublists: int = None

    prefix_execution: str = None
    binary_location: str = None
    pipe_input: str = None
    parallelization: StepExecutionParallelizationParameters = (
        StepExecutionParallelizationParameters()
    )
    failure_policy: StepFailurePolicyParameters = StepFailurePolicyParameters()
    check_backend_availability: bool = False
    resources: StepExecutionResourceParameters = StepExecutionResourceParameters()
    platform: _EPE = _EPE.LOCAL


class StepSettingsArgsParameters(BaseModel):
    flags: List = []
    parameters: Dict = {}


class StepSettingsParameters(BaseModel):
    arguments: StepSettingsArgsParameters = StepSettingsArgsParameters()
    additional: Dict = {}


class StepBase(BaseModel):
    step_id: str
    work_dir: str = None
    type: str = None
    data: StepData = StepData()
    input: StepInputParameters = StepInputParameters()
    writeout: List[StepWriteoutParameters] = []
    execution: StepExecutionParameters = StepExecutionParameters()
    settings: StepSettingsParameters = StepSettingsParameters()

    class Config:
        underscore_attrs_are_private = True

    _logger = PrivateAttr()
    _logger_blank = PrivateAttr()
    _old_wdir = PrivateAttr()
    _workflow_object = PrivateAttr()
    _backend_executor: Executor = PrivateAttr()
    _subtask_container: SubtaskContainer = PrivateAttr()

    def __init__(self, **data):
        super().__init__(**data)

        self._logger_blank = BlankLogger()
        self._old_wdir = os.getcwd()
        self._workflow_object = None
        self._backend_executor = None

        self._logger = StepLogger()
        self._logger_blank = BlankLogger()

    # @staticmethod
    def _make_tmpdir(self):
        if self.work_dir is not None:
            return self.work_dir
        else:
            self.work_dir = tempfile.mkdtemp()
            return self.work_dir

    def _remove_temporary(self, paths):
        if paths is not None:
            if not isinstance(paths, list):
                paths = [paths]
            if (
                self.get_workflow_object() is None
                or self.get_workflow_object().header.global_settings.remove_temporary_files
            ):
                for path in paths:
                    if os.path.isdir(path):
                        shutil.rmtree(path, ignore_errors=True)
                    elif os.path.isfile(path) and os.path.exists(path):
                        os.remove(path)
                    else:
                        self._logger.log(
                            f"Path {path} is neither a valid folder nor file path.",
                            _LE.WARNING,
                        )
            else:
                self._logger.log(
                    f"Keeping {len(paths)} temporary file(s) / folder(s): {', '.join(paths)}",
                    _LE.DEBUG,
                )

    @staticmethod
    def _move_to_temp_dir() -> str:
        cur_tmp_dir = tempfile.mkdtemp()
        os.chdir(cur_tmp_dir)
        return cur_tmp_dir

    @staticmethod
    def _move_to_dir(path: str):
        os.chdir(path)

    def _restore_working_dir(self):
        os.chdir(self._old_wdir)

    def execute(self):
        raise NotImplementedError

    def get_compound_by_name(self, name: str) -> Compound:
        for compound in self.data.compounds:
            if compound.get_name() == name:
                return compound

    def get_compounds(self) -> List[Compound]:
        return self.data.compounds

    def get_generic(self) -> GenericContainer:
        return self.data.generic

    def clone_compounds(self) -> List[Compound]:
        return [deepcopy(comp) for comp in self.data.compounds]

    def process_write_out(self):
        for writeout in self.writeout:
            writeout_handler = WriteOutHandler(config=writeout)
            writeout_handler.set_data(self.data)
            # attach workflow data at this point
            writeout_handler.set_workflow_data(self.get_workflow_object().workflow_data)
            writeout_handler.write()

    def get_compound_stats(self) -> Tuple[int, int, int]:
        n_comp = len(self.get_compounds())
        n_enum = len(unroll_enumerations(self.get_compounds()))
        n_conf = len(unroll_conformers(self.get_compounds()))
        return n_comp, n_enum, n_conf

    def generate_input(self):
        preparator = InputPreparator(
            workflow=self.get_workflow_object(), logger=self._logger
        )
        self.data, self.work_dir = preparator.generate_input(
            step_input=self.input, step_type=self.type
        )

        # check for a perturbation map for fep workflows
        self._logger.log(
            f"Loaded {len(self.data.compounds)} compounds and {len(self.data.generic.get_flattened_files())} generic data fields for step {self.get_step_id()}.",
            _LE.DEBUG,
        )

    def set_workflow_object(self, workflow_object):
        self._workflow_object = workflow_object

    def get_workflow_object(self):
        return self._workflow_object

    def get_perturbation_map(self) -> PerturbationMap:
        return self._workflow_object.workflow_data.perturbation_map

    def get_step_id(self) -> str:
        return self.step_id

    def set_step_id(self, step_id: str):
        self.step_id = step_id

    def _initialize_backend(self, executor: Callable):
        if self.execution.platform == _EPE.SLURM:
            self._backend_executor = SlurmExecutor(
                prefix_execution=self.execution.prefix_execution,
                binary_location=self.execution.binary_location,
                cores=self.execution.resources.cores,
                tasks=self.execution.resources.tasks,
                partition=self.execution.resources.partition,
                time=self.execution.resources.time,
                mem=self.execution.resources.mem,
                modules=self.execution.resources.modules,
                other_args=self.execution.resources.other_args,
                additional_lines=self.execution.resources.additional_lines,
                gres=self.execution.resources.gres,
            )
        else:
            self._backend_executor = executor(
                prefix_execution=self.execution.prefix_execution,
                binary_location=self.execution.binary_location,
            )

    def _unroll_compounds(
        self,
        compounds: List[Compound],
        level: str = _SBE.WRITEOUT_COMP_CATEGORY_CONFORMERS,
    ) -> List[Conformer]:
        # TODO: move this to step_base or merge with methods from compound itself

        all_conformers = []
        for comp in compounds:
            for enum in comp.get_enumerations():
                if level == _SBE.WRITEOUT_COMP_CATEGORY_ENUMERATIONS:
                    all_conformers.append(enum)
                elif level == _SBE.WRITEOUT_COMP_CATEGORY_CONFORMERS:
                    for conf in enum:
                        all_conformers.append(conf)
        return all_conformers

    def write_conformers(self, path: str):
        """Convenience function for frequent conformer coordinate write-out. Better to use the WriteOutHandler class."""
        compounds_copy = self.clone_compounds()
        params = {
            _SBE.WRITEOUT_CONFIG: {
                _SBE.WRITEOUT_COMP: {
                    _SBE.WRITEOUT_COMP_CATEGORY: _SBE.WRITEOUT_COMP_CATEGORY_CONFORMERS
                },
                _SBE.WRITEOUT_DESTINATION: {
                    _SBE.WRITEOUT_DESTINATION_RESOURCE: path,
                    _SBE.WRITEOUT_DESTINATION_TYPE: _SBE.WRITEOUT_DESTINATION_TYPE_FILE,
                    _SBE.WRITEOUT_DESTINATION_FORMAT: _SBE.FORMAT_SDF,
                },
            }
        }
        writeout_handler = WriteOutHandler(**params)
        writeout_handler.set_data(StepData(compounds=compounds_copy))
        writeout_handler.write()

    def write_generic_by_extension(self, path: str, ext: str, join=True):
        """writes all files of a specific file type to the specified directory, retaining original files names"""
        for file in self.data.generic.get_files_by_extension(ext):
            file.write(path, join=join)

    def write_generic_by_name(self, path, name: str):
        file = self.data.generic.get_file_by_name(name)
        file.write(path)

    def _check_backend_availability(self, strict=True):
        if self._backend_executor is None:
            raise Exception(
                "Cannot check backend availability before initialization is complete."
            )

        if self.execution.check_backend_availability:
            if not self._backend_executor.is_available():
                if strict:
                    raise StepFailed(
                        f"Cannot initialize backend for step {self.step_id} - abort."
                    )
                else:
                    self._logger.log(
                        f"Backend availability check failed, proceeding anyways.",
                        _LE.WARNING,
                    )
            else:
                self._logger.log(f"Checked backend availability - valid.", _LE.DEBUG)

    def _input_object_valid(self, obj) -> bool:
        if obj.get_molecule() is None or not isinstance(obj.get_molecule(), Chem.Mol):
            self._logger.log(
                f"Object {obj.get_index_string()} skipped - no valid molecule.",
                _LE.WARNING,
            )
            return False
        return True

    def _input_object_empty(self, obj) -> bool:
        if obj.empty():
            self._logger.log(
                f"Object {obj.get_index_string()} is skipped (empty).", _LE.WARNING
            )
            return True
        return False

    # TODO: REMOVE THIS FUNCTION (see: write_molecule_to_sdf())
    def _prepare_temp_input(self, tmp_dir: str, molecule: Chem.Mol) -> str:
        _, tmp_sdf_path = gen_tmp_file(suffix=".sdf", dir=tmp_dir)
        if molecule is None or not isinstance(molecule, Chem.Mol):
            raise ValueError(
                "Function requires input attribute to be an RDkit molecule."
            )
        writer = Chem.SDWriter(tmp_sdf_path)
        writer.write(molecule)
        writer.close()
        self._logger.log(f"Wrote input molecule to file {tmp_sdf_path}.", _LE.DEBUG)
        return tmp_sdf_path

    def _get_sublists(self, get_first_n_lists: int = None) -> List[List[Subtask]]:
        number_cores = self._get_number_cores()

        # decide how to slice the ligand list depending on whether a maximum length is defined or not
        if self.execution.parallelization.max_length_sublists is not None:
            slice_size = min(
                max(self.execution.parallelization.max_length_sublists, 1),
                len(self._subtask_container),
            )
            return self._subtask_container.get_sublists(
                partitions=None,
                slice_size=slice_size,
                get_first_n_lists=get_first_n_lists,
            )
        else:
            # split the ligands into as many cores as available
            partitions = min(number_cores, len(self._subtask_container))
            return self._subtask_container.get_sublists(
                partitions=partitions,
                slice_size=None,
                get_first_n_lists=get_first_n_lists,
            )

    def _get_number_cores(self):
        # prepare the parallelization and set the number of cores to be used
        cores = self.execution.parallelization.jobs
        if cores == 0:
            cores = 1
        elif cores < 0:
            # subtract the number of cores (neg. value, thus add up) from total number of cores, e.g. -1 will
            # use all available cores minus 1
            cores = multiprocessing.cpu_count() + cores
        return cores

    def _print_log_file(self, path: str):
        if os.path.isfile(path):
            with open(path, "r") as log_file:
                log_file_raw = log_file.readlines()
                self._logger.log(f"Printing log file {path}:\n", _LE.DEBUG)
                for line in log_file_raw:
                    self._logger_blank.log(line.rstrip("\n"), _LE.DEBUG)
                self._logger_blank.log("", _LE.DEBUG)
                self._logger.log("--- End file", _LE.DEBUG)

    def _add_data_to_generic(self, file, data, extension=None):
        """Write data from arbitrary file to generic container class"""
        file_name = file.split("/")[-1]
        # file types where they can be passed as arguments in a subsequent step
        # TODO: this is not maintainable!
        file_tag = (
            True
            if file.endswith((".gro", "topol.top", "tpr", "fmp", "edge"))
            else False
        )
        file = GenericData(
            file_name=file_name, file_data=data, argument=file_tag, extension=extension
        )
        self.data.generic.add_file(file)

    def _parse_output(
        self,
        tmp_dir,
        exclusion_list=(
            "#",
            "AC",
            "AC0",
            "INF",
            "hashed",
            "metadata",
            "timekeys",
            "000000",
        ),
    ):
        """Generic method for parsing generic writeout, can be overwritten in child classes"""
        self.data.generic.clear_file_dict()
        file_list = [os.path.join(tmp_dir, f) for f in os.listdir(tmp_dir)]
        for file in file_list:
            if os.path.isfile(file) and not file.endswith(exclusion_list):
                try:
                    with open(file, "r") as f:
                        data = f.read()
                    binary = False
                except UnicodeDecodeError:
                    with open(file, "rb") as f:
                        data = f.read()
                    binary = True
                # work out if we handle the data or just the path to it on disk
                file_size = os.stat(file).st_size
                if file_size > float(_SBE.FILE_SIZE_THRESHOLD.value):
                    # do not write to the dict - file is too large to store in memory
                    _, tmp_path = gen_tmp_file(suffix="." + str(file).split(".")[-1])
                    self._logger.log(
                        f"Large file detected, storing at {tmp_path}", _LE.INFO
                    )
                    if binary:
                        with open(tmp_path, "wb") as f:
                            f.write(data)
                    else:
                        with open(tmp_path, "w") as f:
                            f.write(data)
                    data = tmp_path

                self._add_data_to_generic(file, data)
                self._logger.log(f"Stored data for file {file}", _LE.DEBUG)
            elif os.path.isdir(file):
                tmp_dir = mkdtemp()
                copy_tree(file, tmp_dir)
                self._add_data_to_generic(file=file, data=tmp_dir, extension="dir")

                # we have picked up a directory, we want the entire contents copied somewhere

    def _wait_until_file_generation(
        self,
        path,
        path_log=None,
        interval_sec=1,
        maximum_sec=None,
        success_strings: set = set(),
        fail_strings: set = set(),
    ) -> bool:
        # TODO: Refactor that without breaking the Glide dependency.
        counter = 0
        while not os.path.exists(path):
            # wait for an interval
            time.sleep(interval_sec)
            counter = counter + 1

            # if a Glide logfile path has been specified, check, whether critical messages indicating an abort are there
            # note, that we return "True" to indicate that the "file generation" has nevertheless been completed
            if path_log is not None:
                if any_in_file(path_log, fail_strings):
                    self._logger.log(
                        f"A critical error occurred in sublist execution.", _LE.WARNING
                    )
                    self._print_log_file(path_log)
                    return True
                if any_in_file(path_log, success_strings):
                    # log file indicates job is done; give a bit of leeway to ensure the writing is done
                    time.sleep(3)
                    break

            # if there's time left, proceed
            if maximum_sec is not None and counter * interval_sec >= maximum_sec:
                break
        if os.path.exists(path):
            return True
        else:
            return False

    def _log_execution_progress(self):
        number_tasks_done = len(self._subtask_container.get_done_tasks())
        number_tasks_total = len(self._subtask_container.subtasks)
        self._logger.log(
            get_progress_bar_string(number_tasks_done, number_tasks_total, length=65),
            _LE.INFO,
        )

    def _log_result(self, result: CompletedProcess):
        """
        logs stdout from completed process to file
        """
        for line in result.stdout.split("\n"):
            self._logger_blank.log(line, _LE.DEBUG)

    def get_additional_setting(self, key: str, default: str = None):
        """
        Query settings.additional with the key, if not set use the default
        """
        return (
            self.settings.additional[key]
            if key in self.settings.additional.keys()
            else default
        )

    def get_topol(self) -> GromacsState:
        if not self.data.generic.get_file_names_by_extension("pkl"):
            return self.data.gmx_state
        else:
            return self.load_topol(
                self.data.generic.get_argument_by_extension("pkl", rtn_file_object=True)
            )
