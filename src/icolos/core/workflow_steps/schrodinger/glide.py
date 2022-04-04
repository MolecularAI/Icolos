import gzip
import os
import shutil
import tempfile
from copy import deepcopy
from typing import List, Tuple

from pydantic import BaseModel
from rdkit import Chem

from icolos.core.workflow_steps.schrodinger.base import StepSchrodingerBase
from icolos.utils.execute_external.glide import GlideExecutor
from icolos.utils.execute_external.schrodinger import SchrodingerExecutor
from icolos.utils.general.files_paths import any_in_file, gen_tmp_file

from icolos.core.containers.compound import Conformer

from icolos.utils.enums.program_parameters import SchrodingerExecutablesEnum, GlideEnum
from icolos.utils.enums.step_enums import StepGlideEnum, StepBaseEnum
from icolos.core.workflow_steps.step import _LE
from icolos.utils.general.parallelization import SubtaskContainer, Parallelizer, Subtask
from icolos.utils.general.strings import stringify


class GlideSupportEnum:

    GLIDE_INPUTBLOCK_COMMASEPARATED = [
        "CONSTRAINT_GROUP"
    ]  # define list of block keys which are to have commas
    GLIDE_INPUTBLOCK_VALUEQUOTED = [
        "FEATURE"
    ]  # define list of block keys, where values are to be put
    # into double quotation marks

    GLIDE_TG_WAIT_INTERVAL = "wait_interval_seconds"
    GLIDE_TG_WAIT_LIMIT = "wait_limit_seconds"

    # try to find the internal value and return
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    # prohibit any attempt to set any values
    def __setattr__(self, key, value):
        raise ValueError("No changes allowed.")


_SBE = StepBaseEnum
_EE = GlideEnum()
_SGE = StepGlideEnum()
_SEE = SchrodingerExecutablesEnum()
_GSE = GlideSupportEnum()


class StepGlide(StepSchrodingerBase, BaseModel):

    _schrodinger_executor: SchrodingerExecutor = None

    class Config:
        underscore_attrs_are_private = True

    def __init__(self, **data):
        super().__init__(**data)

        # initialize the executors and test availability
        self._initialize_backend(executor=GlideExecutor)
        self._check_backend_availability()

        self._schrodinger_executor = SchrodingerExecutor(
            prefix_execution=self.execution.prefix_execution,
            binary_location=self.execution.binary_location,
        )

    def _get_scores_from_conformer(self, conformer: Chem.Mol) -> Tuple[float, float]:
        return (
            float(conformer.GetProp(_SGE.GLIDE_DOCKING_SCORE)),
            float(conformer.GetProp(_SGE.GLIDE_GSCORE)),
        )

    def _set_docking_score(self, conformer: Chem.Mol) -> bool:
        try:
            docking_score, g_score = self._get_scores_from_conformer(conformer)
        except KeyError:
            return False
        conformer.SetProp(_SBE.ANNOTATION_TAG_DOCKING_SCORE, str(docking_score))
        conformer.SetProp(_SBE.ANNOTATION_TAG_G_SCORE, str(g_score))
        return True

    def _generate_temporary_input_output_files(
        self, batch: List[List[Subtask]]
    ) -> Tuple[List[str], List[str], List[str], List[str]]:
        tmp_output_dirs = []
        tmp_input_mae_paths = []
        tmp_output_sdf_paths = []
        tmp_output_maegz_paths = []

        for next_subtask_list in batch:
            # generate temporary input files and output directory
            cur_tmp_output_dir = tempfile.mkdtemp()
            _, cur_tmp_sdf = gen_tmp_file(suffix=".sdf", dir=cur_tmp_output_dir)
            _, cur_tmp_mae = gen_tmp_file(suffix=".mae", dir=cur_tmp_output_dir)

            # write-out the temporary input file
            writer = Chem.SDWriter(cur_tmp_sdf)
            one_written = False
            for subtask in next_subtask_list:
                enumeration = subtask.data
                mol = deepcopy(enumeration.get_molecule())
                if mol is not None:
                    mol.SetProp("_Name", enumeration.get_index_string())
                    one_written = True
                    writer.write(mol)
            writer.close()
            if one_written is False:
                self._remove_temporary(cur_tmp_output_dir)
                continue

            # translate the SDF into a MAE file
            self._translate_SDF_to_MAE(
                sdf_path=cur_tmp_sdf,
                mae_path=cur_tmp_mae,
                executor=self._schrodinger_executor,
            )

            # add the path to which "_dock_subjob()" will write the result SDF
            _, output_sdf_path = gen_tmp_file(
                suffix="_result.sdf", dir=cur_tmp_output_dir
            )
            _, output_maegz_path = gen_tmp_file(
                suffix="_result.maegz", dir=cur_tmp_output_dir, text=False
            )
            tmp_output_sdf_paths.append(output_sdf_path)
            tmp_output_maegz_paths.append(output_maegz_path)
            tmp_input_mae_paths.append(cur_tmp_mae)
            tmp_output_dirs.append(cur_tmp_output_dir)
        return (
            tmp_output_dirs,
            tmp_input_mae_paths,
            tmp_output_sdf_paths,
            tmp_output_maegz_paths,
        )

    def _all_keywords(self) -> dict:
        """Returns joined keywords from JSON and from .in file (if specified)."""

        keywords = {}

        # keywords from maestro file; they can be overwritten by explicitly set values from the "configuration" block
        maestro_in_file = deepcopy(
            self.settings.additional.get(_SGE.MAESTRO_IN_FILE, None)
        )
        if maestro_in_file is not None:
            with open(maestro_in_file[_SGE.MAESTRO_IN_FILE_PATH], "rt") as f:
                keywords_from_file = self._parse_maestro_in_file(f.readlines())
                keywords.update(keywords_from_file)

        # Add keywords from advanced_glide_keywords
        # (they are keywords with file paths),
        # skipping keywords that are None.
        # Also skip maestro file - that's not a keyword.
        # TODO: This is legacy code from DockStream's implementation, which was necessary to accommodate the GUI.
        #       Remove?
        # if self.parameters.advanced_glide_keywords is not None:
        #    adv_kw = stringify({
        #        k: v
        #        for k, v in self.parameters.advanced_glide_keywords.dict().items()
        #        if v is not None and k not in {'maestro_file'}
        #    })
        #    keywords.update(adv_kw)

        # Add "ordinary" keywords, overwriting existing ones.
        json_keywords = stringify(
            deepcopy(self.settings.additional.get(_SGE.CONFIGURATION, {}))
        )
        keywords.update(
            json_keywords
        )  # Overwrites any keywords that are already present.
        return keywords

    def _configuration_Maestro_reformat(self, configuration: dict):
        # rewrite keyword input file in Maestro format
        maestro_indent = "    "
        maestro_spacing = "   "

        element_lines = []
        block_lines = []

        for key in configuration.keys():
            if isinstance(configuration[key], str):
                # keyword holds one dictionary (string) only
                element_lines.append(
                    maestro_spacing.join([key, configuration[key] + "\n"])
                )
            elif isinstance(configuration[key], dict):
                # keyword holds a composite block and has no dictionary (e.g. constraints); note, that these must
                # always be at the end of the file
                block_lines.append("\n" + key + "\n")
                block = configuration[key]
                for key_idx, block_key in enumerate(block.keys()):
                    block_value = block[block_key]

                    # if this is a value in certain blocks, put it into double quotation marks as spaces are present
                    if any([x in key for x in _GSE.GLIDE_INPUTBLOCK_VALUEQUOTED]):
                        block_value = '"' + block_value + '"'
                    line = maestro_indent + maestro_spacing.join(
                        [block_key, block_value]
                    )

                    # add comma to block definition, if there are more lines to come and the block requires it
                    # note, that not all blocks in GLIDE require this; in some cases, the comma is already part of
                    # the line (then skip it!)
                    if any([x in key for x in _GSE.GLIDE_INPUTBLOCK_COMMASEPARATED]):
                        if (key_idx + 1) < len(block) and line[-1] != ",":
                            line = line + ","

                    block_lines.append(line + "\n")
            else:
                raise Exception(
                    f"Cannot handle type {type(configuration[key])} in configuration file specification, only use strings and blocks."
                )

        return element_lines, block_lines

    def _write_configuration_to_file(self, configuration: dict, path: str):
        """Function to generate a keyword input file in Maestro format."""

        # call a function that returns the input keywords in Maestro format
        element_lines, block_lines = self._configuration_Maestro_reformat(
            configuration=configuration
        )

        # arrange the elements and blocks
        if path is None:
            _, path = gen_tmp_file(suffix=".in")
        with open(path, mode="w") as f:
            self._logger.log(f"Writing GLIDE input file {path}:\n", _LE.DEBUG)
            for line in element_lines:
                f.write(line)
                self._logger_blank.log(line.rstrip("\n"), _LE.DEBUG)
            for line in block_lines:
                f.write(line)
                self._logger_blank.log(line.rstrip("\n"), _LE.DEBUG)
            self._logger_blank.log("", _LE.DEBUG)
            self._logger.log("--- End file", _LE.DEBUG)

    def _get_time_limit_per_task(self):
        # for "SP" method, it can be expected to that about 90 s / ligand is required at most
        # use a bit extra
        return int(self.settings.additional.get(_SGE.TIME_LIMIT_PER_TASK, 240))

    def _get_path_tmp_results(
        self, glide_pose_outtype: str, base_path: str
    ) -> Tuple[str, str]:
        if glide_pose_outtype == _EE.GLIDE_POSE_OUTTYPE_LIGANDLIB:
            path_tmp_results = os.path.join(
                os.path.dirname(base_path),
                "".join(
                    [
                        os.path.splitext(os.path.basename(base_path))[0],
                        _SGE.GLIDE_SDF_DEFAULT_EXTENSION,
                    ]
                ),
            )
        elif glide_pose_outtype == _EE.GLIDE_POSE_OUTTYPE_POSEVIEWER:
            path_tmp_results = os.path.join(
                os.path.dirname(base_path),
                "".join(
                    [
                        os.path.splitext(os.path.basename(base_path))[0],
                        _SGE.GLIDE_MAEGZ_DEFAULT_EXTENSION,
                    ]
                ),
            )
        else:
            raise NotImplementedError(
                f"Specified out-type {glide_pose_outtype} for Glide not supported."
            )

        path_tmp_log = os.path.join(
            os.path.dirname(base_path),
            "".join([os.path.splitext(os.path.basename(base_path))[0], _SGE.GLIDE_LOG]),
        )
        return path_tmp_results, path_tmp_log

    def _move_result_files(
        self,
        glide_pose_outtype: str,
        path_tmp_results: str,
        path_sdf_results: str,
        path_maegz_results: str,
    ):
        if glide_pose_outtype == _EE.GLIDE_POSE_OUTTYPE_LIGANDLIB:
            if os.path.isfile(path_tmp_results):
                with gzip.open(path_tmp_results, "rb") as fin:
                    with open(path_sdf_results, "wb") as fout:
                        shutil.copyfileobj(fin, fout)
        elif glide_pose_outtype == _EE.GLIDE_POSE_OUTTYPE_POSEVIEWER:
            # as the output is in MAEGZ format, we need to translate it into an SDF (and move the original file to the
            # expected path)
            self._translate_MAE_to_SDF(
                mae_path=path_tmp_results,
                sdf_path=path_sdf_results,
                executor=self._schrodinger_executor,
            )
            os.rename(path_tmp_results, path_maegz_results)
        else:
            raise NotImplementedError(
                f"Specified out-type {glide_pose_outtype} for Glide not supported."
            )

    def _run_subjob(
        self,
        mae_ligand_path,
        path_sdf_results,
        path_maegz_results,
        tmp_output_dir,
        grid_path,
        sublist,
    ):
        # 1) increase the sublist "tries" and set status to "failed"
        _ = [task.increment_tries() for task in sublist]
        _ = [task.set_status_failed() for task in sublist]

        # 2) change to directory, to be able to use relative paths (to compensate for Schrodinger bug with AWS)
        working_dir = os.getcwd()
        os.chdir(tmp_output_dir)

        # 3) get "keywords" dictionary and overwrite necessary values
        #    add "LIGANDFILE" keyword to list of keywords: full path to "mae" formatted ligands
        configuration = self._all_keywords()
        if configuration is None:
            raise ValueError(
                f"You need to specify at least the gridfile path in the configuration for Glide."
            )
        configuration[_EE.GLIDE_LIGANDFILE] = mae_ligand_path

        # set the path to the grid file for this run
        configuration[_EE.GLIDE_GRIDFILE] = grid_path

        # if not set, set the liand pose outtype to "LIGANDLIB" (SDF output without receptor)
        glide_pose_outtype = configuration.get(
            _EE.GLIDE_POSE_OUTTYPE, _EE.GLIDE_POSE_OUTTYPE_LIGANDLIB
        )
        configuration[_EE.GLIDE_POSE_OUTTYPE] = glide_pose_outtype

        # 4) write the keyword-input file for the "Glide" backend; write-out to temporary file
        _, glide_configuration_path = gen_tmp_file(suffix=".in", dir=tmp_output_dir)
        self._write_configuration_to_file(
            configuration=configuration,
            path=glide_configuration_path,
        )

        # 5) wait / sleep until job is completed
        #    Note, that while Glide has an option "-WAIT", this does not seem to work when getting back
        #    data from AWS (probably it ends before copying back the data properly); stay with this solution for now
        path_tmp_results, path_tmp_log = self._get_path_tmp_results(
            glide_pose_outtype=glide_pose_outtype, base_path=glide_configuration_path
        )

        # 5) execute the "Glide" backend
        arguments = self._prepare_glide_arguments(glide_configuration_path)
        execution_result = self._backend_executor.execute(
            command=_EE.GLIDE,
            arguments=arguments,
            check=True,
            location=os.path.dirname(glide_configuration_path),
        )

        # 6) check return code (anything but '0' is bad) and add "stdout" to log file
        time_exceeded = False
        if execution_result.returncode != 0:
            msg = (
                f"Could not dock with Glide, error message: {execution_result.stdout}."
            )
            self._logger.log(msg, _LE.ERROR)
            self._print_log_file(path_tmp_log)
            raise RuntimeError()
        else:
            if (
                self._wait_until_file_generation(
                    path=path_tmp_results,
                    path_log=path_tmp_log,
                    interval_sec=10,
                    maximum_sec=max(
                        self._get_time_limit_per_task() * len(sublist), 300
                    ),
                    success_strings=_EE.GLIDE_LOG_FINISHED_STRINGS,
                    fail_strings=_EE.GLIDE_LOG_FAIL_STRINGS,
                )
                is False
            ):
                time_exceeded = True
                self._logger.log(
                    f"Sublist docking for output file {path_tmp_results} exceeded time limit or failed, "
                    f"all these ligands are ignored in the final write-out. This could mean that none of "
                    f"them could be docked or a runtime error in Glide occured.",
                    _LE.DEBUG,
                )

        # 6) load the log-file (if generated) and check if all went well
        if (
            any_in_file(path_tmp_log, _EE.GLIDE_LOG_SUCCESS_STRING)
            and time_exceeded is False
        ):
            self._logger.log(
                f"Finished sublist (input: {mae_ligand_path}, output: {path_sdf_results}).",
                _LE.DEBUG,
            )
        else:
            self._print_log_file(path_tmp_log)

        # 7) collect the results; Glide outputs the sdf with a given, semi-hard-coded path; extract the sdf file
        self._move_result_files(
            glide_pose_outtype=glide_pose_outtype,
            path_tmp_results=path_tmp_results,
            path_sdf_results=path_sdf_results,
            path_maegz_results=path_maegz_results,
        )

        # 8) revert back to working directory
        os.chdir(working_dir)

    def _prepare_glide_arguments(self, glide_configuration_path: str) -> List[str]:
        # Note, that the first argument is the path to the configuration input file
        # If the number of cores has been set, overwrite "N_JOBS" and parallelize internally and also note
        # that each subjob requires a license; instead start each with "N_JOBS" = 1
        arguments = [glide_configuration_path]

        # copy parameters and overwrite as necessary
        parameters = deepcopy(self.settings.arguments.parameters)
        parameters[_EE.GLIDE_NJOBS] = 1

        if len(self.settings.arguments.flags) > 0:
            for flag in self.settings.arguments.flags:
                # -WAIT leads to issues at times: The process may not return properly
                # (e.g. because of writing problems) and then gets stuck; workaround with waiting
                # for file completion, so remove it if set
                if flag not in [_EE.GLIDE_WAIT]:
                    arguments.append(str(flag))
        if parameters:
            for key in parameters.keys():
                # remove "-WAIT" if set as a parameter, as this leads to instability issues and ignore empty keys
                if key == _EE.GLIDE_WAIT or key == "":
                    continue
                arguments.append(key)
                if parameters[key] is not None and parameters[key] != "":
                    arguments.append(str(parameters[key]))
        return arguments

    def _execute_glide(self, grid_id: str, grid_path: str):
        # TODO: add individual resubmission for failed subtasks
        # get number of sublists in batch and initialize Parallelizer
        glide_parallelizer = Parallelizer(func=self._run_subjob)

        # continue until everything is successfully done or number of retries have been exceeded
        while self._subtask_container.done() is False:
            next_batch = self._get_sublists(get_first_n_lists=self._get_number_cores())

            # generate paths and initialize molecules (so that if they fail, this can be covered)
            (
                tmp_output_dirs,
                tmp_input_mae_paths,
                tmp_output_sdf_paths,
                tmp_output_maegz_paths,
            ) = self._generate_temporary_input_output_files(next_batch)

            # call "token guard" method (only executed, if block is specified in the configuration), which will wait
            # with the execution if not enough tokens are available at the moment
            self._apply_token_guard()

            # execute the current batch in parallel; hand over lists of parameters (will be handled by Parallelizer)
            # also increment the tries and set the status to "failed" (don't do that inside subprocess, as data is
            # copied, not shared!)
            _ = [sub.increment_tries() for element in next_batch for sub in element]
            _ = [sub.set_status_failed() for element in next_batch for sub in element]
            list_grid_path = [grid_path for _ in tmp_input_mae_paths]
            glide_parallelizer.execute_parallel(
                mae_ligand_path=tmp_input_mae_paths,
                path_sdf_results=tmp_output_sdf_paths,
                path_maegz_results=tmp_output_maegz_paths,
                tmp_output_dir=tmp_output_dirs,
                grid_path=list_grid_path,
                sublist=next_batch,
            )

            # parse the output of that particular batch and remove temporary files
            self._parse_glide_output(
                tmp_output_sdf_paths,
                tmp_output_maegz_paths,
                next_batch,
                grid_id,
                grid_path,
            )

            # clean-up
            self._remove_temporary(tmp_output_dirs)

            # print the progress for this execution
            self._log_execution_progress()

    def _log_execution(self, grid_id: str, number_grids: int):
        number_enumerations = 0
        number_conformers = 0
        for compound in self.get_compounds():
            number_enumerations += len(compound)
            for enumeration in compound:
                number_conformers += len(enumeration)
                if len(enumeration) == 0:
                    self._logger.log(
                        f"Enumeration {enumeration.get_index_string()} has no docked poses attached.",
                        _LE.DEBUG,
                    )
        self._logger.log(
            f"Executed Schrodinger/Glide backend for grid {grid_id} (of {number_grids}), now storing a total of {number_conformers} conformers for {number_enumerations} enumerations in {len(self.get_compounds())} compounds.",
            _LE.INFO,
        )

    def _parse_glide_output(
        self,
        tmp_output_sdf_paths: List[str],
        tmp_output_maegz_paths: List[str],
        batch: List[List[Subtask]],
        grid_id: str,
        grid_path: str,
    ):
        # TODO: refactor that (recombine with ligprep parsing?)
        def _update_subtask(sublist: List[Subtask], enum_identifier: str):
            for task in sublist:
                if task.data.get_index_string() == enum_identifier:
                    task.set_status_success()

        def _add_poseviewer_file(conformer: Conformer, maegz_path: str):
            if os.path.isfile(maegz_path) and os.path.getsize(maegz_path) > 0:
                with open(maegz_path, "rb") as f:
                    conformer.add_extra_data(
                        key=_SGE.GLIDE_POSEVIEWER_FILE_KEY, data=f.read()
                    )

        for i in range(len(tmp_output_sdf_paths)):
            # get input and output paths and check the files are there
            path_sdf_results = tmp_output_sdf_paths[i]
            path_maegz_results = tmp_output_maegz_paths[i]
            cur_sublist = batch[i]

            # this is a protection against the case where empty (file size == 0 bytes) files are generated due to
            # a failure during docking
            if (
                not os.path.isfile(path_sdf_results)
                or os.path.getsize(path_sdf_results) == 0
            ):
                continue

            mol_supplier = Chem.SDMolSupplier(path_sdf_results, removeHs=False)
            for mol in mol_supplier:
                if mol is None:
                    continue
                cur_enumeration_name = str(mol.GetProp("_Name"))

                # add the information on the actual grid used
                mol.SetProp(_SBE.ANNOTATION_GRID_ID, str(grid_id))
                mol.SetProp(_SBE.ANNOTATION_GRID_PATH, str(grid_path))
                mol.SetProp(_SBE.ANNOTATION_GRID_FILENAME, os.path.basename(grid_path))

                # if no docking score is attached (i.e. the molecule is a receptor or so, skip it)
                if self._set_docking_score(mol) is not True:
                    continue

                # add molecule to the appropriate ligand
                for compound in self.get_compounds():
                    for enumeration in compound:
                        if enumeration.get_index_string() == cur_enumeration_name:
                            new_conformer = Conformer(
                                conformer=mol,
                                conformer_id=None,
                                enumeration_object=enumeration,
                            )
                            _add_poseviewer_file(
                                conformer=new_conformer, maegz_path=path_maegz_results
                            )
                            enumeration.add_conformer(new_conformer, auto_update=True)
                            _update_subtask(
                                cur_sublist, enum_identifier=cur_enumeration_name
                            )
                            break

    def _sort_conformers(self):
        # sort the conformers (best to worst) and update their names to contain the conformer id
        # -> <compound>:<enumeration>:<conformer_number>
        for compound in self.get_compounds():
            for enumeration in compound:
                enumeration.sort_conformers(
                    by_tag=_SGE.GLIDE_DOCKING_SCORE, reverse=False
                )

    def execute(self):
        # in order to be able to efficiently execute Glide on the enumeration level, all of them have to be unrolled
        # Note: As they retain their respective Compound object, the attribution later on is simple
        all_enumerations = []
        for compound in self.get_compounds():
            all_enumerations = all_enumerations + compound.get_enumerations()
            for enumeration in compound:
                enumeration.clear_conformers()

        # to allow ensemble docking, loop over all provided grid files and annotate the origin of the conformers
        gridfiles = deepcopy(self.settings.additional.get(_SGE.CONFIGURATION, None))[
            _EE.GLIDE_GRIDFILE
        ]
        if not isinstance(gridfiles, list):
            gridfiles = [gridfiles]

        # set grid ids (generate indices, if not specified)
        grid_ids = self.settings.additional.get(_SBE.GRID_IDS, [])
        if len(grid_ids) != len(gridfiles):
            self._logger.log(
                f"There were {len(grid_ids)} grid_ids specified for {len(gridfiles)}, using indices instead.",
                _LE.DEBUG,
            )
            grid_ids = [str(idx) for idx in range(len(gridfiles))]

        for grid_id, grid_path in zip(grid_ids, gridfiles):
            # split into sublists, according to the settings
            self._subtask_container = SubtaskContainer(
                max_tries=self.execution.failure_policy.n_tries
            )
            self._subtask_container.load_data(all_enumerations)

            # execute Glide
            self._execute_glide(grid_id=grid_id, grid_path=grid_path)

            # do the logging
            self._log_execution(grid_id=grid_id, number_grids=len(gridfiles))

        # sort the conformers loaded to the enumerations
        self._sort_conformers()
