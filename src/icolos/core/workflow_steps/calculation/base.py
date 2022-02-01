import numpy as np
import pandas as pd

from pydantic import BaseModel
from rdkit.Chem import AllChem
from typing import List

from icolos.core.containers.compound import Conformer

from icolos.core.workflow_steps.step import StepBase
from icolos.utils.enums.step_enums import StepRMSFilterEnum

_SRF = StepRMSFilterEnum()


class StepCalculationBase(StepBase, BaseModel):
    def __init__(self, **data):
        super().__init__(**data)

    def _get_rms_method(self):
        # there are two modes for the execution: "best" is better, but sometimes has performance issues
        # for larger molecules
        if self.settings.additional[_SRF.METHOD] == _SRF.METHOD_ALIGNMOL:
            return AllChem.AlignMol
        elif self.settings.additional[_SRF.METHOD] == _SRF.METHOD_BEST:
            return AllChem.GetBestRMS
        else:
            raise ValueError(
                f"RMS mode {self.settings.arguments.parameters[_SRF.METHOD]} not supported (either {_SRF.METHOD_ALIGNMOL} or {_SRF.METHOD_BEST})."
            )

    @staticmethod
    def _get_property_values(conformers: List[Conformer], prop: str) -> List[float]:
        return [float(conf.get_molecule().GetProp(prop)) for conf in conformers]

    @staticmethod
    def _calculate_rms_matrix(
        conformers: List[Conformer], rms_method, decimals=3
    ) -> pd.DataFrame:
        n_conf = len(conformers)
        df_rms = pd.DataFrame(np.nan, index=range(n_conf), columns=range(n_conf))
        np.fill_diagonal(df_rms.values, 0)

        for i in range(n_conf - 1):
            for j in range(i + 1, n_conf):
                df_rms.iloc[i, j] = df_rms.iloc[j, i] = np.round(
                    rms_method(
                        conformers[i].get_molecule(), conformers[j].get_molecule()
                    ),
                    decimals=decimals,
                )
        return df_rms
