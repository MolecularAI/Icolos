import os
import tempfile
from typing import Tuple, List
from copy import deepcopy

from pydantic import BaseModel

from icolos.utils.enums.step_enums import StepTurbomoleEnum
from icolos.utils.execute_external.execute import execution_successful
from icolos.utils.execute_external.openbabel import OpenBabelExecutor
from icolos.utils.execute_external.turbomole import TurbomoleExecutor
from icolos.utils.general.convenience_functions import nested_get

from icolos.utils.general.molecules import get_charge_for_molecule

from icolos.core.containers.compound import Conformer, Enumeration

from icolos.utils.enums.program_parameters import OpenBabelEnum
from icolos.utils.enums.program_parameters import TurbomoleEnum
from icolos.utils.enums.compound_enums import ConformerContainerEnum
from icolos.core.workflow_steps.calculation.base import StepCalculationBase
from icolos.core.workflow_steps.step import _LE
from icolos.loggers.logger_utils import log_multiline_string
from icolos.utils.general.files_paths import _FG, check_file_availability, gen_tmp_file

from icolos.utils.general.parallelization import Parallelizer, SubtaskContainer

_OE = OpenBabelEnum()
_EE = TurbomoleEnum()
_COE = ConformerContainerEnum()
_STE = StepTurbomoleEnum()


class StepTurbomole(StepCalculationBase, BaseModel):
    """Step that executes Turbomole.

    Note, that the execution (especially in conjunction with a subsequent cosmo step) is relatively complex.
    (1) Write the conformer as an SDF file to a temporary directory,
    (2) use obabel to translate it to an XYZ file,
    (3) use t2x to make a coord file out of it (input for turbomole; is updated during geometry optimization),
    (4) execute Turbomole, generating a (i) a final coord file and (ii) a trajectory (if specified),
    (5) use x2t to extract the final "snapshot" as an XYZ file and translate it to SDF and
    (6) update the coordinates and tags in the conformers

    IMPORTANT: Keep the "mol.cosmo" file attached to the conformer as additional data for a possible cosmo step."""

    def __init__(self, **data):
        super().__init__(**data)

        # initialize the executor and test availability
        self._initialize_backend(executor=TurbomoleExecutor)
        # TODO: figure out why "module load turbomole/73 && rdfit" sometimes fails (see also below) and
        #       use strict=True after fix; probably, it has to to with $TURBOTMPDIR (all parallel jobs access the same)
        self._check_backend_availability(strict=False)

    def get_original_conformer(self, conformer: Conformer) -> Conformer:
        for compound in self.get_compounds():
            if conformer.get_enumeration_object()._compound_object.get_compound_number() == compound.get_compound_number():
                for enum in compound.get_enumerations():
                    if (
                        enum._enumeration_id
                        == conformer.get_enumeration_object().get_enumeration_id()
                    ):
                        for conf in enum.get_conformers():
                            if conf._conformer_id == conformer._conformer_id:
                                return conf

    def _prepare_tmp_input_directories(
        self, batch: List
    ) -> Tuple[List, List[str], List[str], List[str], List[str], List[str], List[str]]:
        conformers = []
        tmp_dirs = []
        paths_input_sdf = []
        paths_input_xyz = []
        paths_coord = []
        paths_tm_config = []
        paths_cosmo_config = []
        for sublist in batch:
            for element in sublist:  # there is only one
                conformer = element.data
                conformers.append(conformer)
                # 1) generate all temporary paths
                tmp_dir = tempfile.mkdtemp()
                _, path_input_sdf = gen_tmp_file(
                    prefix="tmp_", suffix=".sdf", dir=tmp_dir
                )
                _, path_input_xyz = gen_tmp_file(
                    prefix="tmp_", suffix=".xyz", dir=tmp_dir
                )
                path_coord = os.path.join(tmp_dir, _EE.COORD)

                # 2) write-out the conformers for an enumeration in a SDF file
                conformer.write(path=path_input_sdf)

                # 3) translate the SDF into an XYZ file (using OpenBabel)
                # Note, that all tags are lost here (but the names are not!)
                obabel_executor = OpenBabelExecutor()
                obabel_executor.execute(
                    command=_OE.OBABEL,
                    arguments=[
                        _OE.OBABEL_INPUTFORMAT_SDF,
                        path_input_sdf,
                        _OE.OBABEL_OUTPUT_FORMAT_XYZ,
                        "".join([_OE.OBABEL_O, path_input_xyz]),
                    ],
                    check=True,
                    location=tmp_dir,
                )

                # 4) translate the XYZ to a TM input file ("coord"); "x2t" writes to stdout
                result = self._backend_executor.execute(
                    command=_EE.TM_X2T, arguments=[path_input_xyz], check=True
                )
                with open(path_coord, "w") as file:
                    file.write(result.stdout)

                # 5) add paths
                tmp_dirs.append(tmp_dir)
                paths_input_sdf.append(path_input_sdf)
                paths_input_xyz.append(path_input_xyz)
                paths_coord.append(path_coord)

                tm_path, cosmo_path = self._get_config_paths(conformer)
                paths_tm_config.append(tm_path)
                paths_cosmo_config.append(cosmo_path)

        return (
            conformers,
            tmp_dirs,
            paths_input_sdf,
            paths_input_xyz,
            paths_coord,
            paths_tm_config,
            paths_cosmo_config,
        )

    def _get_config_paths(self, conformer: Conformer) -> Tuple[str, str]:
        try:
            config_dir = self.settings.additional[_EE.TM_CONFIG_DIR]
            config_basename = self.settings.additional[_EE.TM_CONFIG_BASENAME]
            path_cosmo_config = self.settings.additional[_EE.TM_CONFIG_COSMO]
        except KeyError as e:
            raise KeyError("The dir, basename and cosmo paths need to be set.") from e

        charge = str(
            get_charge_for_molecule(
                molecule=conformer._enumeration_object.get_molecule()
            )
        )

        # the path would look like: /opt/Icolos/turbomole_config/b97-3c-ri-d3-def2-mtzvp-int-nosym-charge-1.tm
        path_tm_config = os.path.join(
            config_dir, "".join([config_basename, charge, _EE.TM_CONFIG_ENDING])
        )
        return path_tm_config, path_cosmo_config

    def _execute_define(self, tmp_dir, path_tm_config: str):
        result = self._backend_executor.execute(
            command=_EE.TM_DEFINE,
            arguments=[" ".join(["<", path_tm_config])],
            check=True,
            location=tmp_dir,
        )

        if not execution_successful(result.stderr, _EE.TM_DEFINE_SUCCESS_STRING):
            self._logger.log(
                f"Execution of {_EE.TM_DEFINE} failed for file {path_tm_config}. Error message:",
                _LE.ERROR,
            )
            log_multiline_string(
                logger=self._logger_blank,
                level=_LE.ERROR,
                multi_line_string=result.stdout,
            )

    def _execute_cosmoprep(self, tmp_dir, path_cosmo_config: str):
        result = self._backend_executor.execute(
            command=_EE.TM_COSMOPREP,
            arguments=[" ".join(["<", path_cosmo_config])],
            check=True,
            location=tmp_dir,
        )

        if not execution_successful(result.stderr, _EE.TM_COSMOPREP_SUCCESS_STRING):
            self._logger.log(
                f"Execution of {_EE.TM_COSMOPREP} failed for file {path_cosmo_config}. Error message:",
                _LE.ERROR,
            )
            log_multiline_string(
                logger=self._logger_blank,
                level=_LE.ERROR,
                multi_line_string=result.stdout,
            )

    def _manipulate_control_script(self, path: str):
        # do the following changes to the "control" script in order to generate FINE Cosmo files
        with open(path, "r") as f:
            control = f.readlines()
        new_control = []
        for line in control:
            if line.rstrip("\n") != _EE.CONTROL_COSMO_OUT:
                new_control.append(line)
            else:
                new_control.append("".join([_EE.CONTROL_COSMO_REPLACE, "\n"]))

                # only add this line in case there is no optimization run going on
                if (
                    nested_get(
                        self.settings.additional,
                        [_STE.EXECUTION_MODE],
                        default=_EE.TM_RIDFT,
                    )
                    == _EE.TM_RIDFT
                ):
                    new_control.append("".join([_EE.CONTROL_COSMO_INSERTION, "\n"]))
        with open(path, "w") as f:
            f.writelines(new_control)

    def _get_arguments(self) -> list:
        arguments = []

        # add flags
        for flag in self.settings.arguments.flags:
            arguments.append(flag)

        # flatten the dictionary into a list for command-line execution
        for key in self.settings.arguments.parameters.keys():
            arguments.append(key)
            arguments.append(self.settings.arguments.parameters[key])
        return arguments

    def _execute_run(self, tmp_dir):
        execution_mode = nested_get(
            self.settings.additional, [_STE.EXECUTION_MODE], default=_EE.TM_RIDFT
        )
        result = self._backend_executor.execute(
            command=execution_mode,
            arguments=self._get_arguments(),
            check=True,
            location=tmp_dir,
        )

        if (
            not execution_successful(result.stderr, _EE.TM_RIDFT_SUCCESS_STRING)
            or result.returncode != 0
        ):
            self._logger.log(
                f"Execution of {execution_mode} failed (return code: {result.returncode}). Error message (stdout & stderr):",
                _LE.DEBUG,
            )
            log_multiline_string(
                logger=self._logger_blank,
                level=_LE.DEBUG,
                multi_line_string=result.stdout,
            )
            log_multiline_string(
                logger=self._logger_blank,
                level=_LE.DEBUG,
                multi_line_string=result.stderr,
            )
        return result.returncode

    def _coord2sdf(self, tmp_dir, path_output_xyz: str, path_output_sdf: str):
        # extract the latest snapshot and write it as an XYZ file
        result = self._backend_executor.execute(
            command=_EE.TM_T2X, arguments=[_EE.TM_T2X_C], check=True, location=tmp_dir
        )

        with open(path_output_xyz, "w") as file:
            file.write(result.stdout)

        # translate it to an SDF
        obabel_executor = OpenBabelExecutor()
        obabel_executor.execute(
            command=_OE.OBABEL,
            arguments=[
                _OE.OBABEL_INPUTFORMAT_XYZ,
                path_output_xyz,
                _OE.OBABEL_OUTPUT_FORMAT_SDF,
                "".join([_OE.OBABEL_O, path_output_sdf]),
            ],
            check=True,
        )

    def _parse_output(self, tmp_dirs: List[str], conformers: List[Conformer]):
        results = []
        # load and attach "mol.cosmo" file
        for tmp_dir, conformer in zip(tmp_dirs, conformers):
            result = _STE.SUCCESS
            cosmo_path = os.path.join(tmp_dir, _EE.TM_OUTPUT_COSMOFILE)
            if check_file_availability(path=cosmo_path) != _FG.NOT_GENERATED:
                with open(cosmo_path, "r") as f:
                    file_content = f.readlines()
                conf = self.get_original_conformer(conformer)
                conf.add_extra_data(key=_COE.EXTRA_DATA_COSMOFILE, data=file_content)
                # conformer.add_extra_data(
                #     key=_COE.EXTRA_DATA_COSMOFILE, data=file_content
                # )

            else:
                self._logger.log(
                    f"Could not load cosmo file for {conformer.get_index_string()}, will remove conformer.",
                    _LE.WARNING,
                )
                self._logger.log(
                    f"File {cosmo_path} could not be loaded for {conformer.get_index_string()}.",
                    _LE.DEBUG,
                )
                result = _STE.FAILED

                # set molecule to None removes the 3D coordinates -> will be deleted in the end
                conformer.set_molecule(None)

            # load and attach "coord" file
            coord_file = os.path.join(tmp_dir, _EE.TM_OUTPUT_COORDFILE)
            coord_file_status = check_file_availability(path=coord_file)
            if coord_file_status == _FG.NOT_GENERATED:
                self._logger.log(
                    f"File {coord_file} could not be loaded for {conformer.get_index_string()}.",
                    _LE.DEBUG,
                )
                result = _STE.FAILED
            elif coord_file_status == _FG.GENERATED_EMPTY:
                self._logger.log(
                    f"File {coord_file} is empty for {conformer.get_index_string()}.",
                    _LE.DEBUG,
                )
                result = _STE.FAILED
            elif coord_file_status == _FG.GENERATED_SUCCESS:
                with open(coord_file, "r") as f:
                    file_content = f.readlines()
                    conf = self.get_original_conformer(conformer)
                    conf.add_extra_data(
                        key=_COE.EXTRA_DATA_COORDFILE, data=file_content
                    )

                execution_mode = nested_get(
                    self.settings.additional,
                    [_STE.EXECUTION_MODE],
                    default=_EE.TM_RIDFT,
                )

                # for RIDFT, only the cosmo file is required as coordinates are not updated (no geometry optimization)
                if execution_mode != _EE.TM_RIDFT:
                    path_output_xyz = os.path.join(tmp_dir, _EE.TM_OUTPUT_FINAL_XYZ)
                    path_output_sdf = os.path.join(tmp_dir, _EE.TM_OUTPUT_FINAL_SDF)
                    self._coord2sdf(tmp_dir, path_output_xyz, path_output_sdf)
                    conf = self.get_original_conformer(conformer)
                    conf.update_coordinates(path=path_output_sdf)
            results.append(result)
        return results

    def _clean_failed_conformers(self, enumeration: Enumeration) -> Tuple[int, int]:
        n_conformers_before = len(enumeration.get_conformers())
        enumeration.clean_failed_conformers()
        n_conformers_after = len(enumeration.get_conformers())
        return n_conformers_before, n_conformers_after

    def _run_conformer(
        self,
        conformer: Conformer,
        tmp_dir: str,
        path_tm_config: str,
        path_cosmo_config: str,
    ) -> None:
        self._execute_define(tmp_dir=tmp_dir, path_tm_config=path_tm_config)
        # execute COSMOprep (update "control")
        self._execute_cosmoprep(tmp_dir=tmp_dir, path_cosmo_config=path_cosmo_config)
        # set a necessary environment variable
        os.environ[_EE.TM_TURBOTMPDIR] = tmp_dir
        # update the "control" file
        self._manipulate_control_script(path=os.path.join(tmp_dir, _EE.CONTROL))
        # all ready; start the execution
        self._execute_run(tmp_dir)

        self._logger.log(
            f"Finished Turbomole execution for conformer {conformer.get_index_string()} in directory {tmp_dir}.",
            _LE.DEBUG,
        )

    def _execute_turbomole_parallel(self):
        parallelizer = Parallelizer(func=self._run_conformer)
        n = 1

        while self._subtask_container.done() is False:

            next_batch = self._get_sublists(
                get_first_n_lists=self._get_number_cores()
            )  # return n lists of length max_sublist_length
            _ = [sub.increment_tries() for element in next_batch for sub in element]
            _ = [sub.set_status_failed() for element in next_batch for sub in element]

            (
                conformers,
                tmp_dirs,
                paths_input_sdf,
                paths_input_xyz,
                paths_coord,
                paths_tm_configs,
                paths_cosmo_configs,
            ) = self._prepare_tmp_input_directories(next_batch)

            self._logger.log(
                f"Executing Turbomole for batch {n} containing {len(tmp_dirs)} conformers",
                _LE.INFO,
            )

            parallelizer.execute_parallel(
                conformer=conformers,
                tmp_dir=tmp_dirs,
                path_tm_config=paths_tm_configs,
                path_cosmo_config=paths_cosmo_configs,
            )

            results = self._parse_output(tmp_dirs, conformers)
            for sublist, result in zip(next_batch, results):
                # TODO: this only works if max length sublist == 1, fine for now as that is all turbomole can handle
                for task in sublist:
                    if result == _STE.SUCCESS:
                        task.set_status_success()
                    else:
                        task.set_status_failed()
            self._remove_temporary(tmp_dirs)
            n += 1

    def execute(self):
        all_conformers = []
        for compound in self.get_compounds():
            for enumeration in compound.get_enumerations():
                if self._input_object_empty(enumeration):
                    continue
                for conformer in enumeration.get_conformers():
                    # for efficient parallelisation, unroll all conformers
                    conf = deepcopy(conformer)
                    all_conformers.append(conf)

        self.execution.parallelization.max_length_sublists = 1
        self._subtask_container = SubtaskContainer(
            max_tries=self.execution.failure_policy.n_tries
        )
        self._subtask_container.load_data(all_conformers)
        self._execute_turbomole_parallel()
