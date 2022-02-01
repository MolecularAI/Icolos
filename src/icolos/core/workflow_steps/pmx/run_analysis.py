from typing import Dict, List
from icolos.core.containers.perturbation_map import Edge
from icolos.core.workflow_steps.pmx.base import StepPMXBase
from pydantic import BaseModel
from icolos.core.workflow_steps.step import _LE
from icolos.utils.enums.program_parameters import StepPMXEnum
from icolos.utils.execute_external.pmx import PMXExecutor
from icolos.utils.general.parallelization import SubtaskContainer
import numpy as np
import glob
import pandas as pd
import os

_PSE = StepPMXEnum()


class StepPMXRunAnalysis(StepPMXBase, BaseModel):
    """
    Analyses map, returns both a summary and a full results dataframe, written to top level of work_dir
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
            run_func=self.run_analysis, step_id="pmx_run_analysis"
        )
        self.analysis_summary(edges)

    def _run_analysis_script(
        self, analysispath, stateApath, stateBpath, bVerbose=False
    ):
        fA = " ".join(glob.glob("{0}/*xvg".format(stateApath)))
        fB = " ".join(glob.glob("{0}/*xvg".format(stateBpath)))
        oA = "{0}/integ0.dat".format(analysispath)
        oB = "{0}/integ1.dat".format(analysispath)
        wplot = "{0}/wplot.png".format(analysispath)
        o = "{0}/results.txt".format(analysispath)
        args = " ".join(self.settings.arguments.flags)

        cmd = "$PMX analyse -fA {0} -fB {1} -o {2} -oA {3} -oB {4} -w {5} -t {6} -b {7}".format(
            fA, fB, o, oA, oB, wplot, 298, 100
        )
        # subprocess complains that the command is too long
        os.system(cmd)

        if bVerbose == True:
            fp = open(o, "r")
            lines = fp.readlines()
            fp.close()
            bPrint = False
            for l in lines:
                if "ANALYSIS" in l:
                    bPrint = True
                if bPrint == True:
                    print(l, end="")

    def _read_neq_results(self, fname):
        fp = open(fname, "r")
        lines = fp.readlines()
        fp.close()
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
        rowName = "{0}_{1}_{2}".format(edge, wp, r)
        self.results_all.loc[rowName, "val"] = res[2]
        self.results_all.loc[rowName, "err_analyt"] = res[3]
        self.results_all.loc[rowName, "err_boot"] = res[4]
        self.results_all.loc[rowName, "framesA"] = res[0]
        self.results_all.loc[rowName, "framesB"] = res[1]

    def _summarize_results(self, edges):
        bootnum = 1000
        for edge in edges:
            for wp in self.therm_cycle_branches:
                dg = []
                erra = []
                errb = []
                distra = []
                distrb = []
                for r in range(1, self.get_perturbation_map().replicas + 1):
                    rowName = "{0}_{1}_{2}".format(edge, wp, r)
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

                rowName = "{0}_{1}".format(edge, wp)
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

            #### also collect self.results_summary
            rowNameWater = "{0}_{1}".format(edge, "ligand")
            rowNameProtein = "{0}_{1}".format(edge, "complex")
            dg = (
                self.results_all.loc[rowNameProtein, "val"]
                - self.results_all.loc[rowNameWater, "val"]
            )
            erra = np.sqrt(
                np.power(self.results_all.loc[rowNameProtein, "err_analyt"], 2.0)
                - np.power(self.results_all.loc[rowNameWater, "err_analyt"], 2.0)
            )
            errb = np.sqrt(
                np.power(self.results_all.loc[rowNameProtein, "err_boot"], 2.0)
                - np.power(self.results_all.loc[rowNameWater, "err_boot"], 2.0)
            )
            rowName = edge
            self.results_summary.loc[rowName, "val"] = dg
            self.results_summary.loc[rowName, "err_analyt"] = erra
            self.results_summary.loc[rowName, "err_boot"] = errb

        # final write to disk
        self.results_summary.to_csv(os.path.join(self.work_dir, "results_summary.csv"))
        self.results_all.to_csv(os.path.join(self.work_dir, "resultsAll.csv"))

    def analysis_summary(self, edges):
        for edge in edges:
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
                    self._fill_resultsAll(res, edge, wp, r)

        # the values have been collected now
        # let's calculate ddGs
        self._summarize_results(edges)

    def run_analysis(self, jobs: List[str], bVerbose=True):
        for idx, edge in enumerate(jobs):

            for r in range(1, self.get_perturbation_map().replicas + 1):

                # ligand
                wp = "ligand"
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
                self._run_analysis_script(
                    analysispath, stateApath, stateBpath, bVerbose=bVerbose
                )
                # protein
                wp = "complex"
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
                self._run_analysis_script(
                    analysispath, stateApath, stateBpath, bVerbose=bVerbose
                )
