from icolos.core.containers.compound import Conformer
from icolos.core.workflow_steps.step import StepBase
from pydantic import BaseModel

from icolos.utils.enums.program_parameters import SchrodingerExecutablesEnum
from icolos.utils.execute_external.schrodinger import SchrodingerExecutor
from icolos.utils.general.parallelization import Parallelizer
import pandas as pd
import os
from icolos.core.workflow_steps.step import _LE


_SEE = SchrodingerExecutablesEnum




class StepProteinInteraction(StepBase, BaseModel):
    def __init__(self, **data):
        super().__init__(**data)

        self._backend_executor = SchrodingerExecutor(
            prefix_execution=self.execution.prefix_execution,
            binary_location=self.execution.binary_location,
        )
        
    def _compute_interactions(self, conf: Conformer):
        tmp_dir = self._make_tmpdir()
        print(tmp_dir)
        conf.write(os.path.join(tmp_dir, "structure.sdf"), format_="pdb")
        command = _SEE.PROTEIN_INTERACTION
        args = self.get_arguments(defaults={"-outfile": "output.csv", "-structure": "structure.pdb"})
        
        
        self._backend_executor.execute(command=command, arguments=args, check=True, location=tmp_dir)
        
        self._parse_output(tmp_dir, conf)
        
        
    def _parse_output(self, path: str, conf: Conformer):
        # parse the output file, attach the boolean to the conformer object
        df = pd.read_csv(os.path.join(path, "output.csv"))
        
        # attach full interaction profile
        conf.add_extra_data("interaction_summary", df)
        
    def _penalize_docking_score(self, conf: Conformer, penalty: float):
        # take the docking score, add a penalty
        docking_score = conf.get_molecule().GetProp("docking_score")
        conf.get_molecule().SetProp("penalized_docking_score", str(float(docking_score) - penalty))

    def execute(self):
        """Runs schrodinger's protein_interaction_analysis script"""
        # requies structure file + group 1/2 identifications

        # unroll conformers
        all_confs = []
        for comp in self.get_compounds():
            for enum in comp.get_enumerations():
                for conf in enum.get_conformers():
                    all_confs.append(conf)
        
        # attach the interaction information to the conformer
        for conf in all_confs:
            self._compute_interactions(conf)
        
        # attach modified docking score if specific interaction is absent 
        penalty = self._get_additional_setting("penalty", default=1.0)
        for conf in all_confs:
            df = conf.get_extra_data()["interaction_summary"]
            # penalize for every interaction that is not met
            for base, interaction in self.settings.additional.items():
                interact_summary = df.loc[df["Residue"] == base]["Specific Interactions"]
                if not f"hb to {interaction}" in interact_summary.values[0]:
                    self._logger.log("Penalizing docking score!", _LE.DEBUG)
                    self._penalize_docking_score(conf, penalty)
        


        
        
        