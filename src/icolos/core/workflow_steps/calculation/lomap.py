import lomap
from icolos.core.containers.perturbation_map import PerturbationMap
from icolos.core.workflow_steps.step import StepBase
from pydantic import BaseModel
import networkx as nx
import os
from icolos.utils.enums.program_parameters import OpenBabelEnum

from icolos.utils.execute_external.openbabel import OpenBabelExecutor
from rdkit import Chem

_OBE = OpenBabelEnum()


class StepLomap(StepBase, BaseModel):
    _openbabel_executor: OpenBabelExecutor = None

    def __init__(self, **data):
        super().__init__(**data)
        self._openbabel_executor = OpenBabelExecutor()

    def execute(self):
        # compute either mcs or radial graph, end up with networkX object which needs to be parsed into a PerturbationMap object
        tmp_dir = self._make_tmpdir()
        existing_dir = os.getcwd()
        os.chdir(tmp_dir)
        self.write_conformers(os.path.join(tmp_dir, "confs.sdf"))
        # separate into individual files with obabel
        args = ["-isdf", "confs.sdf", "-osdf", "-O", "out.sdf", "-m"]
        self._openbabel_executor.execute(
            command=_OBE.OBABEL, arguments=args, check=True, location=tmp_dir
        )
        os.remove(os.path.join(tmp_dir, "confs.sdf"))

        # now rename each file such that it matches the index string
        for f in os.listdir(tmp_dir):
            with Chem.SDMolSupplier(f) as s:
                header = s[0].GetProp("_Name")
            os.rename(f, f"{header}.sdf")

        topology = self._get_additional_setting("topology")
        radial = True if topology == "radial" else False
        # hub_compound_name = self._get_additional_setting("hub_compound")
        # assert hub_compound_name is not None
        # hub_compound = self.get_compound_by_name(hub_compound_name)
        # if hub_compound:
        #     raise NameError(
        #         f"Hub compound {hub_compound_name} was not found in the compound list!"
        #     )
        # write mols out to sdf files

        dbmol = lomap.DBMolecules(
            tmp_dir, self.execution.parallelization.jobs, output=True, radial=radial
        )
        strict, loose = dbmol.build_matrices()
        strict_np = strict.to_numpy_2D_array()
        loose_np = loose.to_numpy_2D_array()

        # generate the nx graph
        nx_graph = dbmol.build_graph()

        p_map = PerturbationMap(
            compounds=self.get_compounds(),
        )
        p_map.generate_from_lomap_output(
            os.path.join(tmp_dir, "out_score_with_connection.txt")
        )
        # revert to old wd
        os.chdir(existing_dir)
        self.get_workflow_object().set_perturbation_map(p_map)
