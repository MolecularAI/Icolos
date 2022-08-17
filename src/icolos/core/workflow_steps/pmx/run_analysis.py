from typing import Dict, List
from icolos.core.containers.perturbation_map import Edge
from icolos.core.workflow_steps.pmx.base import StepPMXBase
from pydantic import BaseModel
from icolos.core.workflow_steps.step import _LE
from icolos.utils.execute_external.pmx import PMXExecutor
from icolos.utils.general.parallelization import SubtaskContainer
import numpy as np
import glob
import pandas as pd
import os


class StepPMXRunAnalysis(StepPMXBase, BaseModel):
    """
    Analyses map, returns both a summary and a full results dataframe, written to top level of work_dir, and attaches properties to the compound
    """

    results_summary: pd.DataFrame = None
    results_all: pd.DataFrame = None

    class Config:
        arbitrary_types_allowed = True

    def __init__(self, **data):
        super().__init__(**data)

        self._initialize_backend(executor=PMXExecutor)
        self.results_summary = pd.DataFrame()
        self.results_all = pd.DataFrame()

    def execute(self):

        edges = [e.get_edge_id() for e in self.get_edges()]
        self.execution.parallelization.max_length_sublists = 1
        self._subtask_container = SubtaskContainer(
            max_tries=self.execution.failure_policy.n_tries
        )
        self._subtask_container.load_data(edges)
        self._execute_pmx_step_parallel(
            run_func=self.run_analysis,
            step_id="pmx_run_analysis",
            result_checker=self._check_result,
        )
        self.analysis_summary(self.get_edges())
        # reattach compounds from perturbation map to step for writeout
        self.data.compounds = self.get_perturbation_map().compounds

    def _run_analysis_script(self, analysispath, stateApath, stateBpath):
        fA = " ".join(glob.glob("{0}/*xvg".format(stateApath)))
        fB = " ".join(glob.glob("{0}/*xvg".format(stateBpath)))
        oA = "integ0.dat"
        oB = "integ1.dat"
        wplot = "wplot.png"
        o = "results.txt"
        # TODO: at the moment we ignore flags from the command line
        # args = " ".join(self.settings.arguments.flags)

        cmd = "$PMX analyse"
        args = [
            "--quiet",
            "-fA",
            fA,
            "-fB",
            fB,
            "-o",
            o,
            "-oA",
            oA,
            "-oB",
            oB,
            "-w",
            wplot,
            "-t",
            298,
            "-b",
            100,
        ]
        self._backend_executor.execute(
            command=cmd, arguments=args, location=analysispath, check=False
        )

    def _read_neq_results(self, fname):
        try:
            with open(fname, "r") as fp:
                lines = fp.readlines()
        except FileNotFoundError:
            return
        out = []
        for l in lines:
            l = l.rstrip()
            foo = l.split()
            if "BAR: dG" in l:
                out.append(float(foo[-2]))
            elif "BAR: Std Err (bootstrap)" in l:
                out.append(float(foo[-2]))
            elif "BAR: Std Err (analytical)" in l:
                out.append(float(foo[-2]))
            elif "0->1" in l:
                out.append(int(foo[-1]))
            elif "1->0" in l:
                out.append(int(foo[-1]))
        return out

    def _fill_resultsAll(self, res, edge, wp, r):
        try:
            rowName = "{0}_{1}_{2}".format(edge, wp, r)
            self.results_all.loc[rowName, "val"] = res[2]
            self.results_all.loc[rowName, "err_analyt"] = res[3]
            self.results_all.loc[rowName, "err_boot"] = res[4]
            self.results_all.loc[rowName, "framesA"] = res[0]
            self.results_all.loc[rowName, "framesB"] = res[1]
        except IndexError:
            self._logger.log(
                f"Index Error encountered whilst parsing results to summary file for job {edge}/{wp}/{r}",
                _LE.WARNING,
            )

    def _summarize_results(self, edges: List[Edge]):
        bootnum = 1000
        for edge in edges:
            try:
                for wp in self.therm_cycle_branches:
                    dg = []
                    erra = []
                    errb = []
                    distra = []
                    distrb = []
                    for r in range(1, self.get_perturbation_map().replicas + 1):
                        rowName = "{0}_{1}_{2}".format(edge.get_edge_id(), wp, r)
                        dg.append(self.results_all.loc[rowName, "val"])
                        erra.append(self.results_all.loc[rowName, "err_analyt"])
                        errb.append(self.results_all.loc[rowName, "err_boot"])
                        distra.append(
                            np.random.normal(
                                self.results_all.loc[rowName, "val"],
                                self.results_all.loc[rowName, "err_analyt"],
                                size=bootnum,
                            )
                        )
                        distrb.append(
                            np.random.normal(
                                self.results_all.loc[rowName, "val"],
                                self.results_all.loc[rowName, "err_boot"],
                                size=bootnum,
                            )
                        )

                    rowName = "{0}_{1}".format(edge.get_edge_id(), wp)
                    distra = np.array(distra).flatten()
                    distrb = np.array(distrb).flatten()

                    if self.get_perturbation_map().replicas == 1:
                        self.results_all.loc[rowName, "val"] = dg[0]
                        self.results_all.loc[rowName, "err_analyt"] = erra[0]
                        self.results_all.loc[rowName, "err_boot"] = errb[0]
                    else:
                        self.results_all.loc[rowName, "val"] = np.mean(dg)
                        self.results_all.loc[rowName, "err_analyt"] = np.sqrt(
                            np.var(distra) / float(self.get_perturbation_map().replicas)
                        )
                        self.results_all.loc[rowName, "err_boot"] = np.sqrt(
                            np.var(distrb) / float(self.get_perturbation_map().replicas)
                        )

                # also collect self.results_summary
                rowNameWater = "{0}_{1}".format(edge.get_edge_id(), "unbound")
                rowNameProtein = "{0}_{1}".format(edge.get_edge_id(), "bound")
                dg = (
                    self.results_all.loc[rowNameProtein, "val"]
                    - self.results_all.loc[rowNameWater, "val"]
                )
                edge.ddG = dg
                try:
                    edge.node_to.get_conformer().get_molecule().SetProp("ddG", str(dg))
                except AttributeError as e:
                    self._logger.log(
                        f"Could not attach score to mol for edge {edge.get_edge_id()}, defaulting to zero",
                        _LE.WARNING,
                    )
                    edge.node_to.get_conformer().get_molecule().SetProp(
                        "ddG", str(0.00)
                    )

                erra = np.sqrt(
                    np.power(self.results_all.loc[rowNameProtein, "err_analyt"], 2.0)
                    + np.power(self.results_all.loc[rowNameWater, "err_analyt"], 2.0)
                )
                edge.ddG_err = erra
                errb = np.sqrt(
                    np.power(self.results_all.loc[rowNameProtein, "err_boot"], 2.0)
                    + np.power(self.results_all.loc[rowNameWater, "err_boot"], 2.0)
                )
                rowName = edge.get_edge_id()

                self.results_summary.loc[rowName, "lig1"] = edge.get_edge_id().split(
                    "_"
                )[0]
                self.results_summary.loc[rowName, "lig2"] = edge.get_edge_id().split(
                    "_"
                )[1]
                self.results_summary.loc[rowName, "val"] = dg
                self.results_summary.loc[rowName, "err_analyt"] = erra
                self.results_summary.loc[rowName, "err_boot"] = errb

            except KeyError as e:
                print(f"Error in generating summary, error was {e}")

    def analysis_summary(self, edges: List[Edge]):
        edge_ids = [e.get_edge_id() for e in edges]
        for edge in edge_ids:
            for r in range(1, self.get_perturbation_map().replicas + 1):
                for wp in self.therm_cycle_branches:
                    analysispath = "{0}/analyse{1}".format(
                        self._get_specific_path(
                            workPath=self.work_dir, edge=edge, wp=wp
                        ),
                        r,
                    )
                    resultsfile = "{0}/results.txt".format(analysispath)
                    res = self._read_neq_results(resultsfile)
                    if res is not None:
                        self._fill_resultsAll(res, edge, wp, r)

        # the values have been collected now
        # let's calculate ddGs
        self._summarize_results(edges)
        # compare with experimental results automatically if provided
        try:
            if "exp_results" in self.settings.additional.keys() and os.path.isfile(
                self.settings.additional["exp_results"]
            ):
                exp_data = pd.read_csv(
                    self.settings.additional["exp_results"],
                    converters={"Ligand": lambda x: str(x).split(".")[0]},
                )
                # compute the experimental ddG and append to resultsSummary
                node_data = self.get_perturbation_map().node_df
                self.results_summary["exp_ddG"] = self.results_summary.apply(
                    lambda x: np.array(
                        self.compute_exp_ddG(
                            x["lig1"], x["lig2"], node_data=node_data, exp_data=exp_data
                        )
                    ),
                    axis=1,
                )
        except Exception as e:
            self._logger.log(
                f"Failed to compute experimental results, error was: {e}", _LE.WARNING
            )
        # final write to disk
        self.results_summary.to_csv(os.path.join(self.work_dir, "resultsSummary.csv"))
        self.results_all.to_csv(os.path.join(self.work_dir, "resultsAll.csv"))

    def compute_exp_ddG(
        self, lig1: str, lig2: str, node_data: pd.DataFrame, exp_data: pd.DataFrame
    ) -> float:
        """
        Compute the ddG between two ligands from experimental data
        """
        lig1_id = (
            node_data.loc[node_data["hash_id"] == lig1]["node_id"]
            .to_list()[0]
            .replace(" ", "")
        )
        lig2_id = (
            node_data.loc[node_data["hash_id"] == lig2]["node_id"]
            .to_list()[0]
            .replace(" ", "")
        )
        lig1_dG = float(
            exp_data.loc[exp_data["Ligand"] == lig1_id]["Exp. ΔG"].tolist()[0]
        )
        lig2_dG = float(
            exp_data.loc[exp_data["Ligand"] == lig2_id]["Exp. ΔG"].tolist()[0]
        )
        return lig2_dG - lig1_dG

    def run_analysis(self, jobs: List[str]):
        for idx, edge in enumerate(jobs):

            for r in range(1, self.get_perturbation_map().replicas + 1):

                # ligand
                wp = "unbound"
                analysispath = "{0}/analyse{1}".format(
                    self._get_specific_path(workPath=self.work_dir, edge=edge, wp=wp),
                    r,
                )
                os.makedirs(analysispath, exist_ok=True)
                stateApath = self._get_specific_path(
                    workPath=self.work_dir,
                    edge=edge,
                    wp=wp,
                    state="stateA",
                    r=r,
                    sim="transitions",
                )
                stateBpath = self._get_specific_path(
                    workPath=self.work_dir,
                    edge=edge,
                    wp=wp,
                    state="stateB",
                    r=r,
                    sim="transitions",
                )
                self._run_analysis_script(analysispath, stateApath, stateBpath)
                # protein
                wp = "bound"
                analysispath = "{0}/analyse{1}".format(
                    self._get_specific_path(workPath=self.work_dir, edge=edge, wp=wp),
                    r,
                )
                os.makedirs(analysispath, exist_ok=True)
                stateApath = self._get_specific_path(
                    workPath=self.work_dir,
                    edge=edge,
                    wp=wp,
                    state="stateA",
                    r=r,
                    sim="transitions",
                )
                stateBpath = self._get_specific_path(
                    workPath=self.work_dir,
                    edge=edge,
                    wp=wp,
                    state="stateB",
                    r=r,
                    sim="transitions",
                )
                self._run_analysis_script(analysispath, stateApath, stateBpath)

    def _check_result(self, batch: List[List[str]]) -> List[List[bool]]:
        """
        Look in each hybridStrTop dir and check the output pdb files exist for the edges
        """
        output_files = ["integ0.dat", "integ1.dat", "results.txt", "wplot.png"]
        analyse_folders = [
            f"analyse{i}" for i in range(1, self.get_perturbation_map().replicas + 1)
        ]
        results = []
        for subjob in batch:
            subjob_results = []
            for job in subjob:
                subjob_results.append(
                    all(
                        [
                            os.path.isfile(
                                os.path.join(self.work_dir, job, branch, folder, f)
                            )
                            for f in output_files
                            for branch in self.therm_cycle_branches
                            for folder in analyse_folders
                        ]
                    )
                )
            results.append(subjob_results)
        return results
