import json
import os
import re
from copy import deepcopy
from tempfile import mkdtemp
from typing import List, Tuple
from pydantic import BaseModel

from icolos.core.containers.compound import Conformer
from icolos.core.workflow_steps.calculation.base import StepCalculationBase
from icolos.utils.enums.program_parameters import JazzyEnum
from icolos.utils.enums.step_enums import StepJazzyEnum
from icolos.core.workflow_steps.step import _LE
from icolos.utils.execute_external.jazzy import JazzyExecutor
from icolos.utils.general.files_paths import gen_tmp_file
from icolos.utils.general.parallelization import SubtaskContainer, Parallelizer

_SJE = StepJazzyEnum()
_JE = JazzyEnum()

_all_jazzy_commands = [
    _JE.VEC,
    _JE.VIS
]

_all_jazzy_properties = [
    _JE.RESULT_DGA,
    _JE.RESULT_DGP,
    _JE.RESULT_DGTOT,
    _JE.RESULT_SA,
    _JE.RESULT_SDC,
    _JE.RESULT_SDX
]


class JazzyAdditional(BaseModel):
    command: str = _JE.VEC


class StepJazzy(StepCalculationBase, BaseModel):

    jazzy_additional: JazzyAdditional = None

    def __init__(self, **data):
        super().__init__(**data)

        self._initialize_backend(executor=JazzyExecutor)
        self._check_backend_availability()

        # initialize the additional settings
        self.jazzy_additional = JazzyAdditional(**self.settings.additional)
        if self.jazzy_additional.command not in _all_jazzy_commands:
            raise ValueError(f"Jazzy command {self.jazzy_additional.command} unknown - abort.")

    def _prepare_batch(self, batch) -> Tuple:
        tmp_dirs = []
        input_smiles = []
        output_files = []
        conformers = []
        for next_subtask_list in batch:
            tmp_dir = mkdtemp()
            tmp_dirs.append(tmp_dir)
            for subtask in next_subtask_list:
                _, tmp_out_path = gen_tmp_file(suffix=".out", dir=tmp_dir)
                output_files.append(tmp_out_path)
                conformer = subtask.data
                conformers.append(conformer)

                # TODO: not really efficient for SMILES, but keep it (?) on conformer level in anticipation of
                #       structural input that differs per conformer
                input_smile = conformer.get_enumeration_object().get_smile()
                input_smiles.append(input_smile)
        return tmp_dirs, input_smiles, output_files, conformers

    def _prepare_arguments(self, settings: List) -> List:

        # add flags
        for flag in self.settings.arguments.flags:
            settings.append(flag)

        # add parameters
        parameters = deepcopy(self.settings.arguments.parameters)

        # flatten the dictionary into a list for command-line execution
        for key in parameters.keys():
            if key in _all_jazzy_commands:
                self._logger.log(
                    f"Use the additional block to specify Jazzy commands, parameter {key} ignored.",
                    _LE.WARNING,
                )
                continue
            settings.append(key)
            settings.append(parameters[key])
        return settings

    def _run_subjob(self, tmp_dir: str, input_smile: str, output_file: str) -> None:
        work_dir = os.getcwd()
        os.chdir(tmp_dir)

        # construct the specified command-line call; only one command can be used at a time
        # e.g. jazzy vec [OPTIONS] SMILES
        arguments = [self.jazzy_additional.command, '"' + input_smile + '"']
        arguments = self._prepare_arguments(arguments)

        result = self._backend_executor.execute(
            command=_JE.JAZZY, arguments=arguments, check=False
        )

        # Jazzy prints the result to stdout -> store it in a temporary file
        with open(output_file, "w") as f:
            f.writelines(result.stdout)

        os.chdir(work_dir)

    def _parse_jazzy_result(
        self, output_files: List[str], conformers: List[Conformer]
    ) -> List:

        results = []

        for output_file, conformer in zip(output_files, conformers):
            # load the JSON string that was captured from stdout and written to the output file
            try:
                with open(output_file) as file:
                    # Jazzy does not output valid JSONs (' instead of "), so we need to replace those
                    # except escaped ones
                    result = file.read().replace("\r", "").replace("\n", "")
                    p = re.compile('(?<!\\\\)\'')
                    result = p.sub('\"', result)
                    result = json.loads(result)
            except:
                self._logger.log(
                    f"Jazzy result for conformer {conformer.get_index_string()} stored in file {output_file} not found - proceeding.",
                    _LE.WARNING
                )
                results.append(_SJE.FAILURE)
                continue

            # attach the properties obtained as tags
            for key in result.keys():
                if key in _all_jazzy_properties:
                    conformer.get_molecule().SetProp(key, str(result[key]))
            results.append(_SJE.SUCCESS)
        return results

    def _execute_jazzy(self):
        jazzy_parallelizer = Parallelizer(func=self._run_subjob)
        n = 1
        while self._subtask_container.done() is False:
            next_batch = self._get_sublists(get_first_n_lists=self._get_number_cores())
            tmp_dirs, input_smiles, output_files, conformers = self._prepare_batch(
                next_batch
            )

            _ = [sub.increment_tries() for element in next_batch for sub in element]
            _ = [sub.set_status_failed() for element in next_batch for sub in element]

            self._logger.log(f"Executing Jazzy for batch {n}.", _LE.DEBUG)

            jazzy_parallelizer.execute_parallel(
                tmp_dir=tmp_dirs, input_smile=input_smiles, output_file=output_files
            )

            results = self._parse_jazzy_result(output_files, conformers)
            for sublist, result in zip(next_batch, results):
                assert len(sublist) == 1
                for task in sublist:
                    if result == _SJE.SUCCESS:
                        task.set_status_success()
                    else:
                        task.set_status_failed()
            n += 1
            self._remove_temporary(tmp_dirs)

    def execute(self):
        all_conformers = []
        for compound in self.get_compounds():
            for enumeration in compound.get_enumerations():
                if enumeration.get_conformers():
                    for conformer in enumeration.get_conformers():
                        all_conformers.append(conformer)
        self.execution.parallelization.max_length_sublists = 1
        self._subtask_container = SubtaskContainer(
            max_tries=self.execution.failure_policy.n_tries
        )
        self._subtask_container.load_data(all_conformers)
        self._execute_jazzy()
        self._logger.log(
            f"Completed execution of Jazzy for {len(all_conformers)} conformers (using their SMILES strings).",
            _LE.INFO,
        )
