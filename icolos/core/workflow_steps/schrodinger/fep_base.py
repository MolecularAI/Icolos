from pydantic import BaseModel
from icolos.core.workflow_steps.schrodinger.base import StepSchrodingerBase
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path
from icolos.utils.enums.step_enums import StepFepPlusEnum
from typing import List
import time
import os
from icolos.core.workflow_steps.step import _LE

_SFE = StepFepPlusEnum()


class StepFEPBase(StepSchrodingerBase, BaseModel):
    """
    Base class containing common functionality for Schrodinger FEP+ workflows
    """

    def __init__(self, **data):
        super().__init__(**data)

    def _parse_output(self, tmp_dir):
        # pick up the final annotated map construction
        self.data.generic.clear_file_dict()
        self._logger.log(f"Reading output map.", _LE.INFO)
        data = None
        counts = 0
        # hold whilst the job data gets written to local fs
        while data is None and counts < 50000:
            try:
                path = [
                    file
                    for file in os.listdir(tmp_dir)
                    if file.endswith(_SFE.FMP_OUTPUT_FILE)
                ]
                assert len(path) == 1
                path = path[0]
                with open(os.path.join(tmp_dir, path), "rb") as f:
                    data = f.read()
            except AssertionError:
                self._logger.log(
                    "Output file has not yet appeared in the file system, sleeping and retrying...",
                    _LE.INFO,
                )
                time.sleep(15)
                counts += 1

        self._add_data_to_generic(path, data)

    def _extract_log_file_data(self, tmp_dir):
        """
        Parses FEP log file to extract edge and node properties
        """
        lines = None
        counts = 0
        # wait whilst job sits in the queue
        while lines is None and counts < 50000:
            try:
                log_file = [
                    file for file in os.listdir(tmp_dir) if file.endswith(_SFE.LOGFILE)
                ]
                assert len(log_file) == 1
                log_file = log_file[0]

                with open(os.path.join(tmp_dir, log_file), "r") as f:
                    lines = f.readlines()

                edge_header_index = [
                    idx for idx, s in enumerate(lines) if _SFE.EDGE_HEADER_LINE in s
                ][-1]
                node_header_index = [
                    idx for idx, s in enumerate(lines) if _SFE.NODE_HEADER_LINE in s
                ][-1]
                end_of_data_index = [
                    idx for idx, s in enumerate(lines) if _SFE.DATA_TERMINUS in s
                ][0]

                edge_data_lines = [
                    line
                    for line in lines[edge_header_index + 3 : node_header_index - 1]
                ]
                node_data_lines = [
                    line
                    for line in lines[node_header_index + 3 : end_of_data_index - 1]
                ]

                self._process_edge_lines(edge_data_lines)
                self._process_node_lines(node_data_lines)

            except AssertionError:
                self._logger.log(
                    "Log file has not yet appeared in the file system, sleeping and retrying...",
                    _LE.INFO,
                )
                time.sleep(15)
                counts += 1

    def _process_node_lines(self, data: List[str]) -> None:
        for entry in data:
            fields = entry.split()
            idx = fields[1]
            dG = fields[2]
            # attach dG tags to compound objects if present
            if self.data.compounds:
                # account for running this step compoundless
                self.data.compounds[int(idx[0])].get_enumerations()[0].get_conformers()[
                    0
                ].get_molecule().SetProp("dG", str(dG))
            self._logger.log(
                f"dG directly from the output file for compound {idx} is {dG} ",
                _LE.INFO,
            )

    def _process_edge_lines(self, edge_data: List[str]) -> None:
        """
        Calibrate dG values using a reference compound and edge ddG from log file output, return dG for each compound
        """

        # caluclate the max ligand index, accounting for ligands that may have been skipped in previous steps, so can't rely on self.get_compounds()
        len_nodes = 0
        for line in edge_data:
            parts = line.split()

            lig_from = int(parts[1].split(":")[0])
            lig_to = int(parts[3].split(":")[0])
            for idx in [lig_from, lig_to]:
                if idx > len_nodes:
                    len_nodes = idx
        len_nodes += 1  # account for zero indexed ligands

        error_matrix = np.zeros((len_nodes, len_nodes))
        ddG_matrix = np.zeros((len_nodes, len_nodes))
        for line in edge_data:
            parts = line.split()
            try:
                # parse the compound info from the log file
                lig_from = int(parts[1].split(":")[0])
                lig_to = int(parts[3].split(":")[0])
                ddG = float(parts[4].split("+-")[0])
                err = float(parts[4].split("+-")[1])
            except ValueError:
                self._logger.log(
                    f"Line: {line} from the logfile contained an unexpected datatype - cannot process this edge - skipping",
                    _LE.WARNING,
                )
                continue

            error_matrix[lig_from, lig_to] = err
            error_matrix[lig_to, lig_from] = err
            ddG_matrix[lig_from, lig_to] = ddG
            ddG_matrix[lig_to, lig_from] = -ddG
        error_matrix = csr_matrix(error_matrix)
        # compute shortest path from one ligand to the anchor
        _, predecessors = shortest_path(
            error_matrix, directed=False, return_predecessors=True, indices=0
        )
        self._construct_dg_per_compound(ddG_matrix, predecessors, error_matrix)

    def _construct_dg_per_compound(
        self, ddG: np.ndarray, predecessors: List, error_matrix: np.ndarray
    ) -> None:
        """
        Calculate the calibrated binding free energy per compound using a reference value
        Attach calcualted dG to compounds
        """
        try:
            ref_dG = self.settings.additional[_SFE.REFERENCE_DG]
        except KeyError:
            self._logger.log(
                "Expected to find a reference dG value for the lead compound, but none was found."
                "Defaulting to 0.00, you will need to apply a manual correction afterwards",
                _LE.WARNING,
            )
            ref_dG = 0.00

        def _calculate_dg(comp_num: int, dG=ref_dG, err=0):
            prev_index = predecessors[comp_num]
            dG += ddG[prev_index, comp_num]
            err += error_matrix[prev_index, comp_num]
            if prev_index != 0:
                _calculate_dg(prev_index, dG=dG, err=err)
            else:
                data = str(round(dG, 2)) + "+-" + str(round(err, 2))
                self.data.compounds[idx].get_enumerations()[0].get_conformers()[
                    0
                ].get_molecule().SetProp("map_dG", data)
                self._logger.log(
                    f"Calculated dG from spanning tree for compound {idx} is {data}",
                    _LE.INFO,
                )

        for comp in self.get_compounds():
            idx = comp.get_compound_number()
            # check whether the compound appeared in the final map
            try:

                if idx == 0:
                    comp.get_enumerations()[0].get_conformers()[
                        0
                    ].get_molecule().SetProp(
                        "map_dG", str(self.settings.additional[_SFE.REFERENCE_DG])
                    )
                if idx != 0:  # skip the reference compound
                    _calculate_dg(idx)
            except IndexError:
                self._logger.log(
                    f"Compound {idx} was not found in the output map, it was likely dropped during the workflow",
                    _LE.WARNING,
                )
                continue
