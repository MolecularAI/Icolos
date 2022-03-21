import os
from copy import deepcopy
from tempfile import mkdtemp
from typing import List, Tuple
from pydantic import BaseModel
from rdkit import Chem

from icolos.core.containers.compound import Conformer
from icolos.core.workflow_steps.calculation.base import StepCalculationBase
from icolos.utils.enums.program_parameters import KallistoEnum, OpenBabelEnum
from icolos.utils.enums.step_enums import StepKallistoEnum
from icolos.core.workflow_steps.step import _LE
from icolos.utils.execute_external.kallisto import KallistoExecutor
from icolos.utils.execute_external.openbabel import OpenBabelExecutor
from icolos.utils.general.files_paths import gen_tmp_file
from icolos.utils.general.icolos_exceptions import StepFailed
from icolos.utils.general.parallelization import SubtaskContainer, Parallelizer

_SKE = StepKallistoEnum()
_KE = KallistoEnum()
_OBE = OpenBabelEnum()

_all_kallisto_commands = [
    _KE.ALP,
    _KE.BONDS,
    _KE.CNS,
    _KE.EEQ,
    _KE.EXS,
    _KE.LIG,
    _KE.PROX,
    _KE.RMS,
    _KE.SORT,
    _KE.STM,
    _KE.VDW,
]


class KallistoAdditional(BaseModel):
    features: List[str] = [_KE.ALP, _KE.BONDS]  # list of features to be obtained


class StepKallisto(StepCalculationBase, BaseModel):

    kallisto_additional: KallistoAdditional = None
    _openbabel_executor: OpenBabelExecutor = None

    def __init__(self, **data):
        super().__init__(**data)

        self._initialize_backend(executor=KallistoExecutor)
        self._check_backend_availability()

        # initialize the executor for all "OpenBabel"
        self._openbabel_executor = OpenBabelExecutor()
        if not self._openbabel_executor.is_available():
            raise StepFailed(
                "Kallisto requires OpenBabel execution, initialization failed."
            )

        # initialize the additional settings
        self.kallisto_additional = KallistoAdditional(**self.settings.additional)

    def _prepare_temp_input(self, tmp_dir: str, molecule: Chem.Mol) -> str:
        tmp_sdf_path = super()._prepare_temp_input(tmp_dir, molecule)

        # Kallisto expects the input either in a mol2 file or an XYZ file
        _, tmp_xyz_path = gen_tmp_file(suffix=".xyz", dir=tmp_dir)

        # translate the
        arguments = [
            tmp_sdf_path,
            _OBE.OBABEL_INPUTFORMAT_SDF,
            _OBE.OBABEL_OUTPUT_FORMAT_XYZ,
            "".join([_OBE.OBABEL_O, tmp_xyz_path]),
        ]
        self._openbabel_executor.execute(
            command=_OBE.OBABEL, arguments=arguments, check=False
        )

        self._logger.log(
            f"Translated input molecule to file {tmp_xyz_path}.", _LE.DEBUG
        )
        return tmp_xyz_path

    def _prepare_batch(self, batch) -> Tuple:
        tmp_dirs = []
        input_files = []
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
                input_xyz_file = self._prepare_temp_input(
                    tmp_dir, conformer.get_molecule()
                )
                input_files.append(input_xyz_file)
        return tmp_dirs, input_files, output_files, conformers

    def _prepare_arguments(self, settings: List) -> List:

        # add flags
        for flag in self.settings.arguments.flags:
            settings.append(flag)

        # add parameters
        parameters = deepcopy(self.settings.arguments.parameters)

        # flatten the dictionary into a list for command-line execution
        for key in parameters.keys():
            if key in _all_kallisto_commands:
                self._logger.log(
                    f"Use the additional block to specify Kallisto commands, parameter {key} ignored.",
                    _LE.WARNING,
                )
                continue
            settings.append(key)
            settings.append(parameters[key])
        return settings

    def _run_subjob(self, tmp_dir: str, input_file: str, output_file: str) -> None:
        work_dir = os.getcwd()
        os.chdir(tmp_dir)

        # construct the specified command-line call
        arguments = []
        for feature in self.kallisto_additional.features:
            if feature not in _all_kallisto_commands:
                self._logger.log(
                    f"Kallisto feature {feature} not supported, will be ignored.",
                    _LE.WARNING,
                )
                continue
            arguments.append(feature)
            arguments.append(input_file)
        arguments = self._prepare_arguments(arguments)

        result = self._backend_executor.execute(
            command=_KE.KALLISTO, arguments=arguments, check=False
        )

        # Kallisto prints the result to stdout -> store it in a temporary file
        with open(output_file, "w") as f:
            f.writelines(result.stdout)

        os.chdir(work_dir)

    def _parse_kallisto_result(
        self, output_files: List[str], conformers: List[Conformer]
    ) -> List:
        def _split_list(inp: List, chunk_size: int) -> List[List[str]]:
            assert len(inp) % chunk_size == 0
            return [
                inp[idx : idx + chunk_size] for idx in range(0, len(inp), chunk_size)
            ]

        results = []
        number_features = len(self.kallisto_additional.features)

        for output_file, conformer in zip(output_files, conformers):
            number_atoms = conformer.get_molecule().GetNumAtoms()

            # load features from output file
            with open(output_file, "r") as f:
                features_lines = f.readlines()
            features_lines = [line.lstrip().rstrip() for line in features_lines]

            # check, that all features could be calculated for all atoms
            expected_lines = number_features * number_atoms
            if expected_lines != len(features_lines):
                self._logger.log(
                    f"Kallisto result for conformer {conformer.get_index_string()} incomplete ({len(features_lines)} lines instead of the expected {expected_lines}), check {output_file} - proceeding.",
                    _LE.WARNING,
                )
                results.append(_SKE.FAILURE)
                continue

            # group the features and add them to the conformers
            sublists = _split_list(features_lines, number_atoms)
            for feature, sublist in zip(self.kallisto_additional.features, sublists):
                conformer.get_molecule().SetProp(feature, "|".join(sublist))
            results.append(_SKE.SUCCESS)
        return results

    def _execute_kallisto(self):
        kallisto_parallelizer = Parallelizer(func=self._run_subjob)
        n = 1
        while self._subtask_container.done() is False:
            next_batch = self._get_sublists(get_first_n_lists=self._get_number_cores())
            tmp_dirs, input_files, output_files, conformers = self._prepare_batch(
                next_batch
            )

            _ = [sub.increment_tries() for element in next_batch for sub in element]
            _ = [sub.set_status_failed() for element in next_batch for sub in element]

            self._logger.log(f"Executing Kallisto for batch {n}.", _LE.DEBUG)

            kallisto_parallelizer.execute_parallel(
                tmp_dir=tmp_dirs, input_file=input_files, output_file=output_files
            )

            results = self._parse_kallisto_result(output_files, conformers)
            for sublist, result in zip(next_batch, results):
                assert len(sublist) == 1
                for task in sublist:
                    if result == _SKE.SUCCESS:
                        task.set_status_success()
                    else:
                        task.set_status_failed()
            n += 1
            self._remove_temporary(tmp_dirs)

    def execute(self):
        assert len(self.kallisto_additional.features) > 0
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
        self._execute_kallisto()
        self._logger.log(
            f"Completed execution of Kallisto for {len(all_conformers)} conformers, attached the following features: [{', '.join(self.kallisto_additional.features)}].",
            _LE.INFO,
        )
