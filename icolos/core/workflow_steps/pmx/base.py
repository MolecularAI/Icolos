from subprocess import CompletedProcess
from pydantic import BaseModel
from icolos.core.containers.perturbation_map import PerturbationMap
from icolos.core.workflow_steps.step import StepBase
from icolos.utils.enums.program_parameters import GromacsEnum
from icolos.utils.enums.step_enums import StepGromacsEnum
from icolos.utils.execute_external.execute import Executor
from icolos.utils.execute_external.pmx import PMXExecutor
import os
from icolos.utils.general.parallelization import Parallelizer
from icolos.core.workflow_steps.step import _LE
import shutil

_GE = GromacsEnum()
_SGE = StepGromacsEnum()


class StepPMXBase(StepBase, BaseModel):

    _antechamber_executor: Executor = None

    def __init__(self, **data):
        super().__init__(**data)

        self._initialize_backend(executor=PMXExecutor)
        self._check_backend_availability()
        self._antechamber_executor = Executor(prefix_execution=_SGE.AMBERTOOLS_LOAD)

    def _parametrise_protein(
        self,
        protein: str = "protein.pdb",
        path: str = "input/protein",
        output="protein.pdb",
    ):
        # run pdb2gmx on the protein
        pdb2gmx_args = [
            "-f",
            os.path.join(self.work_dir, path, protein),
            "-ignh",
            "-water",
            self.settings.additional["water"],
            "-ff",
            self.settings.additional["forcefield"],
            "-o",
            os.path.join(self.work_dir, path, output),
        ]
        self._gromacs_executor.execute(
            command=_GE.PDB2GMX,
            arguments=pdb2gmx_args,
            check=True,
            location=os.path.join(self.work_dir, path),
        )

    def _clean_pdb_structure(self, tmp_dir: str) -> None:
        files = [file for file in os.listdir(tmp_dir) if file.endswith("pdb")]
        for file in files:
            cleaned_lines = []
            with open(os.path.join(tmp_dir, file), "r") as f:
                lines = f.readlines()
            for line in lines:
                if "ATOM" in line or "HETATM" in line:
                    cleaned_lines.append(line)
            with open(os.path.join(tmp_dir, file), "w") as f:
                f.writelines(cleaned_lines)

    def _parametrisation_pipeline(self, tmp_dir, include_top=False, include_gro=False):
        # main pipeline for producing GAFF parameters for a ligand
        arguments_antechamber = [
            "-i",
            "MOL.pdb",
            "-o",
            "MOL.mol2",
            "-fi",
            "pdb",
            "-fo",
            "mol2",
            "-c",
            "gas",
        ]
        self._logger.log(
            f"Running antechamber on structure {tmp_dir.split('/')[-1]}", _LE.DEBUG
        )
        self._antechamber_executor.execute(
            command=_GE.ANTECHAMBER,
            arguments=arguments_antechamber,
            check=True,
            location=tmp_dir,
        )

        arguments_acpype = [
            os.path.join(_GE.ACPYPE_PATH, _GE.ACPYPE_BINARY),
            "-di",
            "MOL.mol2",
            "-c",
            "gas",
        ]
        self._antechamber_executor.execute(
            command=_GE.PYTHON, arguments=arguments_acpype, location=tmp_dir, check=True
        )
        # search the output dir for the itp file
        acpype_dir = [p for p in os.listdir(tmp_dir) if p.endswith(".acpype")][0]
        itp_file = [
            f
            for f in os.listdir(os.path.join(tmp_dir, acpype_dir))
            if f.endswith("GMX.itp")
        ][0]
        shutil.copyfile(
            os.path.join(tmp_dir, acpype_dir, itp_file),
            # standardized name must be enforced here to make argument
            # parsing easier in subsequent pmx steps
            os.path.join(tmp_dir, "MOL.itp"),
        )
        # for abfe calculations we need the ligand_GMX.top + .gro files as well
        if include_top:
            top_file = [
                f
                for f in os.listdir(os.path.join(tmp_dir, acpype_dir))
                if f.endswith("GMX.top")
            ][0]
            shutil.copyfile(
                os.path.join(tmp_dir, acpype_dir, top_file),
                os.path.join(tmp_dir, top_file),
            )
        if include_gro:
            gro_file = [
                f
                for f in os.listdir(os.path.join(tmp_dir, acpype_dir))
                if f.endswith("GMX.gro")
            ][0]
            shutil.copyfile(
                os.path.join(tmp_dir, acpype_dir, gro_file),
                os.path.join(tmp_dir, gro_file),
            )

    def _execute_pmx_step_parallel(self, run_func, step_id: str):
        """
        Instantiates Icolos's parallelizer object,
        runs the step's execute method,
        checks the reutrn codes i.e. will error if an edge fails
        """
        parallelizer = Parallelizer(func=run_func, collect_rtn_codes=True)
        n = 1
        while self._subtask_container.done() is False:

            next_batch = self._get_sublists(
                get_first_n_lists=self._get_number_cores()
            )  # return n lists of length max_sublist_length
            _ = [sub.increment_tries() for element in next_batch for sub in element]
            _ = [sub.set_status_failed() for element in next_batch for sub in element]
            edges = self._prepare_edges(next_batch)
            # to avoid simultaneous processes logging to the same file, pass the
            self._logger.log(
                f"Executing {step_id} for batch {n}, containing {len(edges)} * {len(edges[0])} edges",
                _LE.INFO,
            )

            rtn_codes = parallelizer.execute_parallel(
                edges=edges,
            )
            assert len(rtn_codes) == len(next_batch)
            for idx, sublist in enumerate(next_batch):
                for task in sublist:  # one edge per sublist
                    if rtn_codes[idx] == 0:
                        task.set_status_success()
                    else:
                        task.set_status_failed()

            n += 1

    def get_arguments(self, defaults: dict = None) -> list:
        """
        Construct pmx-specific arguments from the step defaults,
        overridden by arguments specified in the config file
        """
        arguments = []

        # add flags
        for flag in self.settings.arguments.flags:
            arguments.append(flag)

        # flatten the dictionary into a list for command-line execution
        for key in self.settings.arguments.parameters.keys():
            arguments.append(key)
            arguments.append(self.settings.arguments.parameters[key])

        # add defaults, if not already present
        if defaults is not None:
            for key, value in defaults.items():
                if key not in arguments:
                    arguments.append(key)
                    arguments.append(value)
        return arguments

    def get_edges(self):
        """
        Inspect the map object  passed to the step and extract the edge info
        """

        return self.get_workflow_object().workflow_data.perturbation_map.edges

    def get_nodes(self):
        """
        return the nodes attached to the perturbation map
        """
        return self.get_workflow_object().workflow_data.perturbation_map.nodes

    def _get_line_idx(self, data: list, id_str: str) -> int:
        line = [e for e in data if id_str in e]
        assert len(line) == 1
        line = line[0]
        return data.index(line)

    def _construct_perturbation_map(self, work_dir: str, replicas: int):
        # construct the perturbation map and load in the log file
        log_file = self.data.generic.get_argument_by_extension(
            "log", rtn_file_object=True
        )
        log_file.write(work_dir)
        perturbation_map = PerturbationMap(
            compounds=self.data.compounds,
            protein=self.data.generic.get_argument_by_extension(
                "pdb", rtn_file_object=True
            ),
            replicas=replicas,
        )
        perturbation_map.parse_map_file(
            os.path.join(self.work_dir, log_file.get_file_name())
        )

        self._logger.log(
            f"Initialised perturbation map with {len(perturbation_map.get_nodes())} nodes and {len(perturbation_map.get_edges())} edges",
            _LE.INFO,
        )
        self.get_workflow_object().set_perturbation_map(perturbation_map)

    def _get_line_idx(self, data, id_str) -> int:
        # utility to extract the line index with a specific id string
        line = [e for e in data if id_str in e]
        assert len(line) == 1
        line = line[0]
        return data.index(line)

    def _prepare_edges(self, batch):
        edges = []

        for task in batch:
            task_edges = []
            for element in task:  # for now, only a single element
                task_edges.append(element.data)
            edges.append(task_edges)
        return edges

    def _log_result(self, result: CompletedProcess):
        for line in result.stderr.split("\n"):
            self._logger_blank.log(line, _LE.DEBUG)
