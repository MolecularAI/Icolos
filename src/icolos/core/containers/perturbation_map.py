from typing import Dict, List, Optional
from IPython.lib.display import IFrame
import pandas as pd
from icolos.core.containers.compound import Compound, Conformer
from pyvis.network import Network
from icolos.core.containers.generic import GenericData
from icolos.utils.enums.parallelization import ParallelizationEnum

from icolos.utils.enums.step_enums import StepFepPlusEnum
import os
from pydantic import BaseModel


_SFE = StepFepPlusEnum()
_PE = ParallelizationEnum


class Node(BaseModel):
    """
    Container class for the nodes, wrapper class around a compound object
    """

    class Config:
        arbitrary_types_allowed = True

    node_id: str = None
    node_hash: str = None
    conformer: Conformer = Conformer()
    node_connectivity: List = []

    def __init__(self, **data) -> None:
        super().__init__(**data)

    def get_node_id(self) -> str:
        return self.node_id

    def get_node_color(self):
        # TODO: Expand this so we have different colours for each connectivity number [1,10]
        # this is just a placeholder for now
        thresholds = {i: "c0affe" for i in range(10)}

        num_connections = len(self.node_connectivity)
        return thresholds[num_connections]

    def set_node_id(self, node_id: str):
        self.node_id = node_id

    def get_conformer(self) -> Conformer:
        return self.conformer

    def set_conformer(self, conformer: Conformer) -> None:
        self.conformer = conformer

    def get_node_hash(self) -> str:
        return self.node_hash

    # TODO: add methods here to access connectivity and color attributes


class Edge(BaseModel):
    """
    Simple container class the the edges in the perturbation map, specified entirely by the connected nodes
    """

    class Config:
        arbitrary_types_allowed = True

    node_from: Node = Node()
    node_to: Node = Node()
    total: str = None
    mcs: str = None
    chg: str = None
    softbond: str = None
    min_no_atoms: str = None
    snapCoreRmsd: str = None
    bidirSnapCoreRmsd: str = None
    status: _PE = _PE.STATUS_SUCCESS

    def __init__(self, **data):
        super().__init__(**data)

    def get_source_node_name(self):
        return self.node_from.get_node_hash()

    def get_destination_node_name(self):
        return self.node_to.get_node_hash()

    def get_edge_id(self) -> str:
        # construct the edge ID from the node hashes, separated by '_'
        return f"{self.node_from.get_node_hash()}_{self.node_to.get_node_hash()}"

    def _set_status(self, status: str):
        assert status in [_PE.STATUS_SUCCESS, _PE.STATUS_FAILED]
        self.status = status


class PerturbationMap(BaseModel):
    """Hold a map construction parsed from a csv (probabably from a parsed schrodinger log
    file or something) and provide some utility methods for doing pmx calculations on the edges"""

    class Config:
        arbitrary_types_allowed = True

    nodes: List[Node] = []
    edges: List[Edge] = []
    hash_map: Dict = {}
    compounds: List[Compound] = []
    protein: GenericData = None
    vmap_output: IFrame = None
    replicas: int = 3
    node_df: pd.DataFrame = None
    # prune subsequent edge calculations on error
    strict_execution: str = False

    def __init__(self, **data) -> None:
        super().__init__(**data)

    def _get_line_idx(self, data, id_str) -> int:
        line = [e for e in data if id_str in e]
        assert len(line) == 1
        line = line[0]
        return data.index(line)

    def _get_conformer_by_id(self, comp_id: str) -> Optional[Conformer]:
        # get the compund object based on the ID in the ligand table (compound names). At this stage in the workflow we have only one conformer per enumeration
        try:
            # standard icolos naming conventino
            parts = comp_id.split(":")
            compound_id = parts[0]
            enumeration_id = parts[1]
        except:
            # a non-standard compound name has been used
            compound_id = comp_id
        for compound in self.compounds:
            if compound.get_name().split(":")[0] == compound_id:
                rtn_compound = compound
                enums = rtn_compound.get_enumerations()

                if len(enums) == 1:
                    # easy case, there is only one enumeration, return it's single conformer

                    # at this stage, the docking poses must have been filtered to a single entry
                    # per enumeration (an enumerations should have been filtered on charge state etc.)
                    return enums[0].get_conformers()[0]
                else:
                    # multiple enumerations, must be using Icolos naming or we cannot infer which
                    # enumeration should be used
                    enum = rtn_compound.find_enumeration(
                        enumeration_id=int(enumeration_id)
                    )
                    return enum.get_conformers()[0]

    def parse_map_file(self, file_path) -> None:
        # we need to do some format enforcement here (schrodinger or otherwise)

        with open(file_path, "r") as f:
            data = f.readlines()

        start_edge = self._get_line_idx(data, _SFE.EDGE_HEADER_LINE)
        start_node = self._get_line_idx(data, _SFE.NODE_HEADER_LINE)
        stop_node = self._get_line_idx(data, _SFE.SIMULATION_PROTOCOL)
        edge_info_start = self._get_line_idx(data, _SFE.SIMILARITY)

        # TODO: refactor that
        # clean up the data from schrodinger
        split_data = []
        for line in data:
            split_line = line.split("  ")
            stripped_line = []
            for element in split_line:
                if not element.isspace() and element:
                    stripped_line.append(element.strip())
            split_data.append(stripped_line)

        data = split_data

        map_info = pd.DataFrame(
            data[start_edge + 3 : start_node - 1],
            index=None,
            columns=[
                "Short ID",
                "ligand1 -> ligand2",
                "Bennett ddG",
                "Cycle Closure ddG",
                "Complex dG",
                "Solvent dG",
            ],
        )
        self.node_df = pd.DataFrame(
            data[start_node + 3 : stop_node - 1],
            index=None,
            columns=[
                "hash_id",
                "node_id",
                "Predicted dG",
                "Experimental dG",
                "Predicted Solvation dG",
                "Experimental Solvation dG",
            ],
        )
        edge_info = pd.DataFrame(
            data[edge_info_start + 3 : -1],
            columns=[
                "Short ID",
                "Total",
                "Mcs",
                "Charge",
                "SoftBond",
                "MinimumNumberOfAtom",
                "SnapCoreRmsd",
                "BidirectionSnapCore",
            ],
        ).dropna()
        for hash_id, node_id in zip(self.node_df["hash_id"], self.node_df["node_id"]):
            # map the hashes to the compound IDs
            self.hash_map[hash_id] = node_id
            node = Node(
                node_id=node_id,
                node_hash=hash_id,
                conformer=self._get_conformer_by_id(node_id),
            )
            # generate the Node object to wrap the compound
            self.nodes.append(node)

        for _, edge in edge_info.iterrows():
            edge = Edge(
                node_from=self._get_node_by_hash_id(edge[0].split("_")[0]),
                node_to=self._get_node_by_hash_id(edge[0].split("_")[1]),
                total=edge[1],
                mcs=edge[2],
                chg=edge[3],
                softbond=edge[4],
                min_no_atoms=edge[5],
                snapCoreRmsd=edge[6],
                bidirSnapCoreRmsd=edge[7],
            )
            self.edges.append(edge)
        # process the node info, generate the hash map
        for node in self.nodes:
            self._attach_node_connectivity(node)

    def _attach_node_connectivity(self, node: Node):
        # looks through the constructed edges, returns ids of any edges that have the specified node as one component
        connected_edges = []
        for edge in self.edges:
            if (
                edge.node_from.get_node_hash() == node.node_hash
                or edge.node_to.get_node_hash() == node.node_hash
            ):
                connected_edges.append(edge.get_edge_id())
        node.node_connectivity = connected_edges

    def _get_node_by_node_id(self, node_id: str) -> Node:
        for node in self.nodes:
            if node.node_id == node_id:
                return node

    def _get_node_by_hash_id(self, hash_id: str) -> Node:
        for node in self.nodes:
            if node.node_hash == hash_id:
                return node

    def get_edges(self, alive_only=True) -> List[Edge]:
        if alive_only:
            return [e for e in self.edges if e.status == _PE.STATUS_SUCCESS]
        else:
            return self.edges

    def get_nodes(self) -> List[Node]:
        return self.nodes

    def visualise_perturbation_map(self, write_out_path: str) -> None:
        """method for visualising the data as map with pyvis - Network"""
        vmap = Network(directed=True)
        vmap.barnes_hut()

        # this is not an iterable
        for edge in self.edges:
            vmap.add_node(
                edge.get_source_node_name(), color=edge.node_from.get_node_color()
            )
            vmap.add_node(
                edge.get_destination_node_name(), color=edge.node_to.get_node_color()
            )
            vmap.add_edge(
                source=edge.get_source_node_name(),
                to=edge.get_destination_node_name(),
                length=edge.total,
                label="total: " + str(edge.total),
                title="Mcs: " + str(edge.mcs) + ", SnapCoreRMSD: ",
            )
        self.vmap_output = vmap.show(os.path.join(write_out_path, "vmap.html"))
        # return self.vmap_output

    def get_protein(self) -> GenericData:
        return self.protein

    def get_edge_by_id(self, idx: str) -> Optional[Edge]:
        """
        Performs a look up to return the edge.
        If the idx is a path to a batch file (e.g. for run_simulations), edge is extracted from the filepath
        """
        # handle case where the task is actually a path to a batch script
        if not isinstance(idx, Edge):
            parts = idx.split("/")

            for part in parts:
                for e in self.edges:
                    if part == e:
                        idx = e

        match = [e for e in self.edges if e.get_edge_id() == idx]
        if not match:
            return
        else:
            return match[0]

    def __repr__(self) -> str:
        return f"Icolos Perturbation Map object containing {len(self.edges)} edges and {len(self.nodes)} nodes"

    def __str__(self) -> str:
        return self.__repr__()
