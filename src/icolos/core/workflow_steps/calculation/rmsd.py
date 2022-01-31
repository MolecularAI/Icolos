from typing import List
from pydantic import BaseModel

from icolos.core.containers.compound import Conformer, unroll_conformers
from icolos.utils.enums.step_enums import StepRMSDEnum, StepDataManipulationEnum
from icolos.core.workflow_steps.step import _LE
from icolos.core.workflow_steps.calculation.base import StepCalculationBase

_SR = StepRMSDEnum()
_SDM = StepDataManipulationEnum()


class StepRMSD(StepCalculationBase, BaseModel):
    def __init__(self, **data):
        super().__init__(**data)

        # extend parameters
        if _SR.METHOD not in self.settings.additional.keys():
            self.settings.additional[_SR.METHOD] = _SR.METHOD_ALIGNMOL

    def _calculate_RMSD(self, conformers: List[Conformer]):
        for conf in conformers:
            rmsd_matrix = self._calculate_rms_matrix(
                conformers=[conf] + conf.get_extra_data()[_SDM.KEY_MATCHED],
                rms_method=self._get_rms_method(),
            )

            # use the specified tag name if it is the first value and append an index in case there are more
            for idx, col in enumerate(rmsd_matrix.columns[1:]):
                combined_tag = "".join([_SR.RMSD_TAG, "" if idx == 0 else str(idx)])
                rmsd_value = rmsd_matrix.iloc[[0]][col][0]
                conf.get_molecule().SetProp(combined_tag, str(rmsd_value))
                conf.get_extra_data()[_SDM.KEY_MATCHED][idx].get_molecule().SetProp(
                    combined_tag, str(rmsd_value)
                )

    def execute(self):
        # this assumes that the conformers that are to be matched for the calculation of the RMSD matrix, are attached
        # as a list in a generic data field with a specified key
        conformers = unroll_conformers(compounds=self.get_compounds())
        self._calculate_RMSD(conformers=conformers)
        self._logger.log(
            f"Annotated {len(conformers)} conformers with RMSD values (tag: {_SR.RMSD_TAG}).",
            _LE.INFO,
        )

        # TODO: add a nice pandas DF with the RMSD values to a generic data field
