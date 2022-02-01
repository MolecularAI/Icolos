from rdkit.Chem import Mol
from rdkit.Chem.rdMolDescriptors import CalcNumRings, CalcNumAromaticRings
from pydantic import BaseModel

from icolos.utils.enums.program_parameters import FeatureCounterEnum
from icolos.utils.enums.step_enums import StepFeatureCounterEnum
from icolos.core.workflow_steps.step import _LE
from icolos.core.workflow_steps.calculation.base import StepCalculationBase

_FC = FeatureCounterEnum()
_SFC = StepFeatureCounterEnum()


class StepFeatureCounter(StepCalculationBase, BaseModel):
    def __init__(self, **data):
        super().__init__(**data)

        # extend parameters with defaults
        if _SFC.LEVEL not in self.settings.additional.keys():
            self.settings.additional[_SFC.LEVEL] = _SFC.LEVEL_CONFORMER
            self._logger.log(
                f'No operational level for feature counting specified, defaulting to "{_SFC.LEVEL_CONFORMER}".',
                _LE.INFO,
            )

    def _count_rings(self, mol: Mol):
        number_rings = CalcNumRings(mol)
        mol.SetProp(_FC.PROPERTY_NUM_RINGS, str(number_rings))

    def _count_aromatic_rings(self, mol: Mol):
        number_rings = CalcNumAromaticRings(mol)
        mol.SetProp(_FC.PROPERTY_NUM_AROMATIC_RINGS, str(number_rings))

    def _get_feature_method(self, feature: str):
        if feature == _FC.PROPERTY_NUM_RINGS:
            return self._count_rings
        elif feature == _FC.PROPERTY_NUM_AROMATIC_RINGS:
            return self._count_aromatic_rings
        else:
            raise ValueError(f'Feature "{feature}" not yet supported.')

    def execute(self):
        feature = self.settings.additional[_SFC.FEATURE].lower()
        feature_method = self._get_feature_method(feature=feature)
        level = self.settings.additional[_SFC.LEVEL]
        mol_count = 0
        for compound in self.get_compounds():
            for enumeration in compound.get_enumerations():
                if level == _SFC.LEVEL_ENUMERATION:
                    mol = enumeration.get_molecule()
                    if mol is not None:
                        feature_method(mol)
                        mol_count = mol_count + 1
                elif level == _SFC.LEVEL_CONFORMER:
                    for conformer in enumeration.get_conformers():
                        mol = conformer.get_molecule()
                        if mol is not None:
                            feature_method(mol)
                            mol_count = mol_count + 1
                else:
                    raise ValueError(f'Level "{level}" not supported.')
        self._logger.log(
            f'Counted feature "{feature}" for {mol_count} molecules.', _LE.INFO
        )
