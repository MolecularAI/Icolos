from subprocess import CompletedProcess
from typing import Callable, Dict, List
from pydantic import BaseModel
from icolos.core.containers.compound import Compound, Conformer
from icolos.core.containers.perturbation_map import Node, PerturbationMap
from rdkit.Chem import rdmolops
from icolos.core.containers.compound import Compound, Conformer
from icolos.core.containers.perturbation_map import Node, PerturbationMap
from icolos.core.workflow_steps.step import StepBase
from icolos.utils.enums.parallelization import ParallelizationEnum
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
_PE = ParallelizationEnum


class StepPMXBase(StepBase, BaseModel):

    _antechamber_executor: Executor = None
    _gromacs_executor: Executor = None
    sim_types: List = None
    states: List = None
    therm_cycle_branches: List = None
    run_type: str = None
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

        self._antechamber_executor = Executor()
        self._gromacs_executor = GromacsExecutor(
            prefix_execution=self.execution.prefix_execution
        )
        self.sim_types = ["em", "nvt", "eq", "transitions"]
        self.states = ["stateA", "stateB"]
        # for a normal pmx run this would be "water" and "protein"
        # here we rename for compatibility across abfe/rbfe simulations
        self.therm_cycle_branches = ["ligand", "complex"]

        # simulation setup
        self.run_type = self.get_additional_setting(_SPE.RUN_TYPE, "rbfe")
        self.ff = "amber99sb-star-ildn-mut.ff"
        self.boxshape = self.get_additional_setting(_SPE.BOXSHAPE, "dodecahedron")
        self.boxd = self.get_additional_setting(_SPE.BOXD, 1.5)
        self.water = self.get_additional_setting(_SPE.WATER, "tip3p")
        self.conc = self.get_additional_setting(_SPE.CONC, 0.15)
        self.pname = self.get_additional_setting(_SPE.PNAME, "NaJ")
        self.nname = self.get_additional_setting(_SPE.NNAME, "ClJ")
        self.mdp_prefixes = {
            "em": "em",
            "nvt": "nvt",
            "npt": "npt",
            "eq": "eq",
            "transitions": "ti",
        }

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
        """
        Utility function for getting the right paths from a pmx-type directory structure.  Works for both rbfe and abfe runs
        """
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
        self, simpath, toppath, state, sim_type, empath=None, framestart=0, framestop=1
    ) -> CompletedProcess:
        mdp_path = os.path.join(self.work_dir, "input/mdp")
        mdp_prefix = self.mdp_prefixes[sim_type]

        # TODO: is this a liability? would we ever have more than a single topol file?
        top = "{0}/*.top".format(toppath)
        tpr = "{0}/tpr.tpr".format(simpath)
        mdout = "{0}/mdout.mdp".format(simpath)
        # mdp
        if state == "stateA":
            mdp = "{0}/{1}_l0.mdp".format(mdp_path, mdp_prefix)
        else:
            mdp = "{0}/{1}_l1.mdp".format(mdp_path, mdp_prefix)
        # TODO: deal with nvt/npt for abfe
        # str
        if not sim_type == "transitions":
            if sim_type == "em":
                if self.run_type == "rbfe":
                    inStr = f"{toppath}/ions.pdb"
                elif self.run_type == "abfe":
                    inStr = f"{toppath}/genion.gro"
            elif sim_type in ("eq", "nvt", "npt"):
                inStr = "{0}/confout.gro".format(empath)

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
        elif sim_type == "transitions":
            # significant overhead running 81 different subprocesses, limit to a single call with a very long string (might have to use relative paths)
            grompp_full_cmd = []
            for frame in range(framestart, framestop):
                inStr = "{0}/frame{1}.gro".format(simpath, frame)
                tpr = "{0}/ti{1}.tpr".format(simpath, frame)

                grompp_args = [
                    "gmx grompp",
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
                    "4",
                    "-po",
                    mdout,
                    "&&",
                ]

                grompp_full_cmd += grompp_args
            grompp_full_cmd = " ".join(grompp_full_cmd[:-1])
            result = self._gromacs_executor.execute(
                command=grompp_full_cmd, arguments=[], check=True
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

    def _parametrisation_pipeline(
        self, tmp_dir, conf: Conformer, include_top=False, include_gro=False
    ):
        # main pipeline for producing GAFF parameters for a ligand
        charge_method = self.get_additional_setting(
            key=_SGE.CHARGE_METHOD, default="bcc"
        )
        formal_charge = (
            rdmolops.GetFormalCharge(conf.get_molecule()) if conf is not None else 0
        )
        arguments_acpype = [
            "-di",
            "MOL.pdb",
            "-c",
            charge_method,
            "-a",
            "gaff2",
            "-o",
            "gmx",
            "-n",
            formal_charge,
        ]
        self._logger.log("Generating ligand parameters...", _LE.DEBUG)
        self._antechamber_executor.execute(
            command=_GE.ACPYPE_BINARY,
            arguments=arguments_acpype,
            location=tmp_dir,
            check=True,
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

    def _execute_pmx_step_parallel(
        self,
        run_func: Callable,
        step_id: str,
        result_checker: Callable = None,
        **kwargs,
    ):
        """
        Instantiates Icolos's parallelizer object,
        runs the step's execute method,
        passes any kwargs straight to the run_func
        If result_checker is provided,


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

            if result_checker is not None:
                self._logger.log("Checking execution results...", _LE.DEBUG)
                batch_results = result_checker(jobs)
                good_results = 0
                for task, result in zip(next_batch, batch_results):
                    # returns boolean arrays: False => failed job
                    for subtask, sub_result in zip(task, result):
                        if sub_result == False:
                            subtask.set_status_failed()
                            self._logger.log(
                                f"Warning: job {subtask} failed!", _LE.WARNING
                            )
                            if (
                                self.get_perturbation_map() is not None
                                and self.get_perturbation_map().strict_execution
                                and isinstance(subtask.data, str)
                            ):
                                edge = self.get_perturbation_map().get_edge_by_id(
                                    subtask.data
                                )
                                if edge is not None:
                                    edge._set_status(_PE.STATUS_FAILED)

                        else:
                            subtask.set_status_success()
                            good_results += 1
                self._logger.log(
                    f"EXECUTION SUMMARY: Completed {good_results} jobs successfully (out of {len(jobs) * len(jobs[0])} jobs for step {step_id}",
                    _LE.INFO,
                )

            else:
                self._logger.log(
                    f"Step {step_id} is running without a result checker - failed executions will be ignored!",
                    _LE.WARNING,
                )
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

    def _clean_protein(self):
        existing_itp_files = [
            f
            for f in os.listdir(os.path.join(self.work_dir, "input/protein"))
            if f.endswith("itp") and "Protein" in f
        ]
        if (
            not existing_itp_files
        ):  # no protein itp files, we have a single chain that needs extacting from the top file
            with open(os.path.join(self.work_dir, "input/protein/topol.top"), "r") as f:
                top_lines = f.readlines()

            moltype_line = self._get_line_idx(top_lines, _GE.MOLECULETYPES)

            end_itp_line = self._get_line_idx(top_lines, "; Include water topology")

            moltype = top_lines[moltype_line + 2].split()[0]
            cleaned_top = (
                top_lines[:moltype_line]
                + [f'#include "topol_{moltype}.itp']
                + top_lines[end_itp_line:]
            )

            itp_lines = top_lines[moltype_line:end_itp_line]

            with open(os.path.join(self.work_dir, "input/protein/topol.top"), "w") as f:
                f.writelines(cleaned_top)

            with open(
                os.path.join(self.work_dir, f"input/protein/topol_{moltype}.itp"), "w"
            ) as f:
                f.writelines(itp_lines)

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
            strict_execution=self.get_additional_setting(_SPE.STRICT, default=True),
        )
        perturbation_map.parse_map_file(
            os.path.join(self.work_dir, log_file.get_file_name())
        )

        self._logger.log(
            f"Initialised perturbation map with {len(perturbation_map.get_nodes())} nodes and {len(perturbation_map.get_edges())} edges",
            _LE.INFO,
        )
        self.get_workflow_object().set_perturbation_map(perturbation_map)

    def _prepare_edges(self, batch) -> List[str]:
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

    def _separate_atomtypes(self, lig_path: str) -> None:
        with open(os.path.join(lig_path, "MOL.itp"), "r") as f:
            itp_lines = f.readlines()

        start_idx = self._get_line_idx(itp_lines, _GE.ATOMTYPES)
        stop_index = self._get_line_idx(itp_lines, _GE.MOLECULETYPES)

        atomtype_lines = itp_lines[start_idx:stop_index]
        cleaned_itp_lines = itp_lines[stop_index:]
        with open(os.path.join(lig_path, "MOL.itp"), "w") as f:
            f.writelines(cleaned_itp_lines)

        # process the atomtype lines to remove the bondtype
        # col causes gmx to complain
        cleaned_atomtype_lines = []
        for line in atomtype_lines:
            parts = line.split()
            if len(parts) > 5:
                cleaned_parts = [parts[0]] + parts[2:] + ["\n"]
                cleaned_atomtype_lines.append(" ".join(cleaned_parts))
        with open(os.path.join(lig_path, "ffMOL.itp"), "w") as f:
            f.writelines(cleaned_atomtype_lines)

    def _parametrise_nodes(self, jobs):
        if isinstance(jobs, list):
            node = jobs[0]
        else:
            node = jobs
        if isinstance(node, Node):
            node_id = node.get_node_hash()
            conf = node.conformer
        elif isinstance(node, Compound):
            # in abfe we pass compounds here not edges
            node_id = node.get_index_string()
            conf = node.get_enumerations()[0].get_conformers()[0]
        else:
            raise NotImplementedError(f"Cannot parametrize object of type {type(node)}")
        lig_path = os.path.join(self.work_dir, "input", "ligands", node_id)
        os.makedirs(lig_path, exist_ok=True)
        conf.write(os.path.join(lig_path, "MOL.sdf"), format_="pdb")

        # clean the written pdb, remove anything except hetatm/atom lines
        self._clean_pdb_structure(lig_path)
        # now run ACPYPE on the ligand to produce the topology file
        self._parametrisation_pipeline(lig_path, conf=conf)

        # produces MOL.itp, need to separate the atomtypes directive out into ffMOL.itp for pmx
        # to generate the forcefield later
        self._separate_atomtypes(lig_path)
