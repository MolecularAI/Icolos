import multiprocessing
import os
import tempfile
from typing import List
from pydantic import BaseModel
from icolos.core.workflow_steps.step import StepBase
from sklearn.ensemble import RandomForestRegressor
import numpy as np
from icolos.core.workflow_steps.step import _LE
from icolos.utils.enums.step_enums import StepBaseEnum
from icolos.utils.general.convenience_functions import nested_get
from icolos.utils.enums.step_initialization_enum import StepInitializationEnum
from rdkit import Chem
from dscribe.descriptors import SOAP
from dscribe.kernels import REMatchKernel
from sklearn.preprocessing import normalize
from ase import io

_IE = StepInitializationEnum()


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

    def compute_soap_vectors(self, mol: Chem.Mol) -> np.ndarray:
        tmp_dir = tempfile.mkdtemp()
        Chem.rdmolfiles.MolToXYZFile(mol, os.path.join(tmp_dir, "mol.xyz"))

        atoms = io.read(os.path.join(tmp_dir, "mol.XYZ"))

        soap_desc = SOAP(
            species=["C", "H", "O", "N", "F", "Cl"],
            rcut=5,
            nmax=8,
            lmax=6,
            crossover=True,
        )

        return soap_desc.create(atoms)
