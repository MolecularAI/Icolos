from subprocess import CompletedProcess
from typing import Dict, List
from pydantic import BaseModel
from icolos.core.containers.perturbation_map import PerturbationMap
from icolos.core.workflow_steps.step import StepBase
from icolos.utils.enums.program_parameters import GromacsEnum, StepPMXEnum
from icolos.utils.enums.step_enums import StepGromacsEnum
from icolos.utils.execute_external.execute import Executor
from icolos.utils.execute_external.gromacs import GromacsExecutor
import os
from icolos.utils.general.parallelization import Parallelizer
from icolos.core.workflow_steps.step import _LE
import shutil
import glob

_GE = GromacsEnum()
_SGE = StepGromacsEnum()
_SPE = StepPMXEnum()


class StepPMXBase(StepBase, BaseModel):

    _antechamber_executor: Executor = None
    _gromacs_executor: Executor = None
    sim_types: List = None
    states: List = None
    therm_cycle_branches: List = None
    ff: str = None
    boxshape: str = None
    boxd: float = None
    water: str = None
    conc: float = None
    pname: str = None
    nname: str = None
    mdp_prefixes: Dict = None

    def __init__(self, **data):
        super().__init__(**data)

        self._antechamber_executor = Executor(prefix_execution=_SGE.AMBERTOOLS_LOAD)
        self._gromacs_executor = GromacsExecutor(
            prefix_execution=self.execution.prefix_execution
        )
        self.sim_types = ["em", "eq", "transitions"]
        self.states = ["stateA", "stateB"]
        self.therm_cycle_branches = ["water", "protein"]

        # simulation setup, these should not change
        self.ff = "amber99sb-star-ildn-mut.ff"
        self.boxshape = self.get_setting(_SPE.BOXSHAPE, "dodecahedron")
        self.boxd = self.get_setting(_SPE.BOXD, 1.5)
        self.water = self.get_setting(_SPE.WATER, "tip3p")
        self.conc = self.get_setting(_SPE.CONC, 0.15)
        self.pname = self.get_setting(_SPE.PNAME, "NaJ")
        self.nname = self.get_setting(_SPE.NNAME, "ClJ")
        self.mdp_prefixes = {"em": "em", "eq": "eq", "transitions": "ti"}

    def get_setting(self, key: str, default: str):
        """
        Query settings.additional with the key, if not set use the default
        """
        return (
            self.settings.additional[key]
            if key in self.settings.additional.keys()
            else default
        )

    def _get_specific_path(
        self,
        workPath=None,
        edge=None,
        bHybridStrTop=False,
        wp=None,
        state=None,
        r=None,
        sim=None,
    ):
        if edge == None:
            return workPath
        edgepath = "{0}/{1}".format(workPath, edge)

        if bHybridStrTop == True:
            hybridStrPath = "{0}/hybridStrTop".format(edgepath)
            return hybridStrPath

        if wp == None:
            return edgepath
        wppath = "{0}/{1}".format(edgepath, wp)

        if state == None:
            return wppath
        statepath = "{0}/{1}".format(wppath, state)

        if r == None:
            return statepath
        runpath = "{0}/run{1}".format(statepath, r)

        if sim == None:
            return runpath
        simpath = "{0}/{1}".format(runpath, sim)
        return simpath

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

    def _prepare_single_tpr(
        self, simpath, toppath, state, sim_type, empath=None, frameNum=0
    ) -> CompletedProcess:
        mdp_path = os.path.join(self.work_dir, "input/mdp")
        mdp_prefix = self.mdp_prefixes[sim_type]

        top = "{0}/topol.top".format(toppath)
        tpr = "{0}/tpr.tpr".format(simpath)
        mdout = "{0}/mdout.mdp".format(simpath)
        # mdp
        if state == "stateA":
            mdp = "{0}/{1}_l0.mdp".format(mdp_path, mdp_prefix)
        else:
            mdp = "{0}/{1}_l1.mdp".format(mdp_path, mdp_prefix)
        # str
        if sim_type == "em":
            inStr = "{0}/ions.pdb".format(toppath)
        elif sim_type == "eq":
            inStr = "{0}/confout.gro".format(empath)
        elif sim_type == "transitions":
            inStr = "{0}/frame{1}.gro".format(simpath, frameNum)
            tpr = "{0}/ti{1}.tpr".format(simpath, frameNum)

        grompp_args = [
            "-f",
            mdp,
            "-c",
            inStr,
            "-r",
            inStr,
            "-p",
            top,
            "-o",
            tpr,
            "-maxwarn",
            4,
            "-po",
            mdout,
        ]
        result = self._gromacs_executor.execute(
            command=_GE.GROMPP, arguments=grompp_args, check=True
        )

        self._clean_backup_files(simpath)
        return result

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

    def _execute_pmx_step_parallel(self, run_func, step_id: str, **kwargs):
        """
        Instantiates Icolos's parallelizer object,
        runs the step's execute method,
        passes any kwargs straight to the run_func


        """
        parallelizer = Parallelizer(func=run_func)
        n = 1
        while self._subtask_container.done() is False:

            next_batch = self._get_sublists(
                get_first_n_lists=self._get_number_cores()
            )  # return n lists of length max_sublist_length
            _ = [sub.increment_tries() for element in next_batch for sub in element]
            _ = [sub.set_status_failed() for element in next_batch for sub in element]
            jobs = self._prepare_edges(next_batch)
            self._logger.log(
                f"Executing {step_id} for batch {n}, containing {len(jobs)} * {len(jobs[0])} jobs",
                _LE.INFO,
            )

            parallelizer.execute_parallel(jobs=jobs, **kwargs)

            # TODO: find a reliable way to sort this, ideally by inspecting log files
            for element in next_batch:
                for subtask in element:
                    subtask.set_status_success()

            self._log_execution_progress()
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

    def _prepare_edges(self, batch):
        edges = []

        for task in batch:
            task_edges = []
            for element in task:
                task_edges.append(element.data)
            edges.append(task_edges)
        return edges

    def _log_result(self, result: CompletedProcess):
        for line in result.stderr.split("\n"):
            self._logger_blank.log(line, _LE.DEBUG)

    def _clean_backup_files(self, path):
        toclean = glob.glob("{0}/*#".format(path))
        for clean in toclean:
            os.remove(clean)
