import os
import tempfile
from typing import List
from pydantic import BaseModel
import pandas as pd
from icolos.core.containers.compound import Compound, Enumeration
from icolos.core.workflow_steps.step import StepBase
from sklearn.ensemble import RandomForestRegressor
import numpy as np
from icolos.core.workflow_steps.step import _LE
from icolos.utils.enums.step_enums import StepActiveLearningEnum, StepBaseEnum
from icolos.utils.general.convenience_functions import nested_get
from icolos.utils.enums.step_initialization_enum import StepInitializationEnum
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
from rdkit.Chem import Mol

_IE = StepInitializationEnum()
_SALE = StepActiveLearningEnum()


class ActiveLearningBase(StepBase, BaseModel):
    def __init__(self, **data):
        super().__init__(**data)

    def greedy_acquisition(
        self,
        estimator: RandomForestRegressor,
        X: np.ndarray,
        previous_idx,
        n_instances: int,
    ) -> np.ndarray:
        """
        Implement greedy acquisition strategy, return the n_samples best scores
        """
        try:
            predictions = estimator.predict(X)
        except:
            self._logger.log(
                "Estimator is not fitted, defaulting to random predictions", _LE.INFO
            )
            #     # if not initialized, generate random docking scores (absolute)
            predictions = np.random.uniform(12, 0, len(X))

        # zero those predictions we've seen before
        for idx in previous_idx:
            predictions[idx] = 0
        # smaller before n_instances, largest after
        sorted_preds = np.argpartition(predictions, -n_instances, axis=0)[-n_instances:]
        return sorted_preds

    def _initialize_oracle_step_from_dict(self, step_conf: dict) -> StepBase:
        # note this is a bit of a hack to get around a circular import, we can't use the main util
        _STE = StepBaseEnum
        step_type = nested_get(step_conf, _STE.STEP_TYPE, default=None)
        step_type = None if step_type is None else step_type.upper()
        if step_type in _IE.STEP_INIT_DICT.keys():
            return _IE.STEP_INIT_DICT[step_type](**step_conf)
        else:
            raise ValueError(
                f"Backend for step {nested_get(step_conf, _STE.STEPID, '')} unknown."
            )

    def construct_fingerprints(self, library: pd.DataFrame):
        # add morgan FPs
        library[_SALE.MORGAN_FP] = library.apply(
            lambda x: np.array(
                GetMorganFingerprintAsBitVect(x[_SALE.MOLECULE], 2, nBits=2048),
                dtype=np.float32,
            ),
            axis=1,
        )

        library[_SALE.IDX] = [i for i in range(len(library))]

        # library[_SALE.SOAP_VECTOR] = library.apply(
        #     lambda x: get_soap_vector(x[_SALE.MOLECULE], max_length=max_atom_count),
        #     axis=1,
        # )
        # construct pytorch_geomtetric Graph object based on some hand-computed descriptors
        # library[_SALE.GRAPH] = library.apply(
        #     lambda x: create_graph(
        #         x[_SALE.MOLECULE],
        #         x[self.get_additional_setting(_SALE.CRITERIA)],
        #     ),
        #     axis=1,
        # )

        return library

    def _initialize_oracle(self, compound_list: List[pd.Series]) -> List[StepBase]:
        # list of step configs
        base_oracle_config = self.settings.additional["oracle_config"]
        oracle_steps = []
        for step in base_oracle_config:
            oracle_steps.append(self._initialize_oracle_step_from_dict(step))

        # manually attach the compound objects to the oracle's lead step
        # subsequent steps should take their input from the the first step.
        for idx, compound in enumerate(compound_list):
            cmp = Compound(compound_number=idx)
            cmp.add_enumeration(
                Enumeration(
                    compound_object=cmp,
                    smile=compound[_SALE.SMILES],
                    molecule=compound[_SALE.MOLECULE],
                )
            )
            oracle_steps[0].data.compounds.append(cmp)

        return oracle_steps

    def query_oracle(self, compound_list: List[Mol]) -> List[Compound]:
        """
        Interface function with the oracle method, in the most likely case this is ligprep + docking

        Takes the requested compounds and runs them through the oracle workflow, returns the final compounds with annotations

        Notes:
        This could be an arbitrarily complex workflow, but the only thing that's going to change is the compounds.
        """
        oracle_steps = self._initialize_oracle(compound_list)
        # we have a fully initialized step with the compounds loaded.  Execute them
        for idx, step in enumerate(oracle_steps):
            # for subsequent steps we will need to read in from the previous one
            if idx != 0:
                step.generate_input()
            step.execute()
            step.process_write_out()

        # retrieve compounds from the final step
        final_compounds = oracle_steps[-1].data.compounds
        return final_compounds

    def _extract_final_scores(
        self, compounds: List[Compound], criteria: str, highest_is_best: bool = False
    ) -> np.ndarray:
        """
        Takes a list of compound objects from the oracle and extracts the best score based on the provided criteria
        """
        top_scores = []
        for comp in compounds:
            scores = []
            for enum in comp.get_enumerations():
                for conf in enum.get_conformers():
                    scores.append(float(conf._conformer.GetProp(criteria)))

            # if docking generated no conformers
            if not scores:
                scores.append(0.0)

            best_score = max(scores) if highest_is_best else min(scores)
            top_scores.append(best_score)

        return np.absolute(top_scores)
    
    def check_additional(self, key, val=True) -> bool:

        if (
            key in self.settings.additional.keys()
            and self.settings.additional[key] == val
        ):
            return True
        return False
