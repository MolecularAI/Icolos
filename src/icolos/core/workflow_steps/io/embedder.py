from copy import deepcopy

from pydantic import BaseModel
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

from icolos.core.containers.compound import Conformer
from icolos.utils.general.icolos_exceptions import StepFailed
from icolos.utils.enums.step_enums import StepEmbeddingEnum
from icolos.core.workflow_steps.io.base import StepIOBase

from icolos.core.workflow_steps.step import _LE
from icolos.utils.general.convenience_functions import *
from icolos.utils.smiles import to_mol

_SEE = StepEmbeddingEnum()


class StepEmbedding(StepIOBase, BaseModel):
    def __init__(self, **data):
        super().__init__(**data)

        # extend parameters with defaults
        if _SEE.EMBED_AS not in self.settings.additional.keys():
            self.settings.additional[_SEE.EMBED_AS] = _SEE.EMBED_AS_ENUMERATIONS
            self._logger.log(
                f'No embedding level specified, defaulting to "{_SEE.EMBED_AS_ENUMERATIONS}".',
                _LE.INFO,
            )

    def _smile_to_molecule(self, smile: str) -> Chem.Mol:
        mol = to_mol(smile)
        if mol is None:
            self._logger.log(
                f"The smile {smile} could not be transformed into a molecule and will be skipped.",
                _LE.WARNING,
            )
        return mol

    def _embed_with_RDKit(self, smile: str, parameters: dict) -> Chem.Mol:
        molecule = self._smile_to_molecule(smile)

        # deactivate logger to suppress "missing Hs messages"
        RDLogger.DisableLog("rdApp.*")
        try:
            embed_code = AllChem.EmbedMolecule(
                molecule, randomSeed=42, useRandomCoords=True
            )
        except:
            self._logger.log(
                f'Could not embed molecule with SMILES "{smile}", critical error in "RDkit".',
                _LE.WARNING,
            )
            return None

        status = 0
        if embed_code != -1:
            status = AllChem.UFFOptimizeMolecule(molecule, maxIters=600)
            if status == 1:
                self._logger.log(
                    f"The 3D coordinate generation of molecule {smile} did not converge in time.",
                    _LE.WARNING,
                )
        else:
            self._logger.log(
                f"Could not embed molecule {smile} - no 3D coordinates have been generated.",
                _LE.WARNING,
            )
        RDLogger.EnableLog("rdApp.*")

        # add hydrogens to the molecule (if specified)
        if nested_get(parameters, [_SEE.RDKIT_PROTONATE], default=True):
            molecule = Chem.AddHs(molecule, addCoords=True)

        if embed_code != -1 and status == 0:
            return molecule
        else:
            return None

    def _get_embedding_method(self, parameters: dict) -> str:
        method = nested_get(parameters, [_SEE.METHOD], default=None)
        if method is None:
            error = "Embedding method not set."
            self._logger.log(error, _LE.ERROR)
            raise StepFailed(error)
        return method.upper()

    def _embed_molecule(self, smile: str, parameters: dict) -> Chem.Mol:
        method = self._get_embedding_method(parameters)
        if method == _SEE.METHOD_RDKIT:
            return self._embed_with_RDKit(smile, parameters)
        else:
            self._logger.log(
                f"Specified embedding method {method} not available.", _LE.ERROR
            )

    def execute(self):
        # TODO: REFACTOR
        parameters = deepcopy(self.settings.arguments.parameters)
        embed_as = self.settings.additional[_SEE.EMBED_AS]
        for compound in self.get_compounds():
            if embed_as == _SEE.EMBED_AS_ENUMERATIONS:
                enum_buffer = deepcopy(compound.get_enumerations())
                compound.clear_enumerations()
                for enumeration in enum_buffer:
                    enumeration.clear_molecule()
                    enumeration.clear_conformers()
                    molecule = self._embed_molecule(
                        smile=enumeration.get_smile(), parameters=parameters
                    )
                    if molecule is not None:
                        enumeration.set_molecule(molecule)
                        compound.add_enumeration(enumeration)
                self._logger.log(
                    f"Embedding for compound {compound.get_index_string()} (name: {compound.get_name()}) completed ({len(compound)} of {len(enum_buffer)} enumerations successful).",
                    _LE.INFO,
                )
            elif embed_as == _SEE.EMBED_AS_CONFORMERS:
                # TODO: double-check this bit
                for enumeration in compound.get_enumerations():
                    enumeration.clear_conformers()
                    molecule = self._embed_molecule(
                        smile=enumeration.get_smile(), parameters=parameters
                    )
                    if molecule is not None:
                        conformer = Conformer(
                            conformer=molecule, enumeration_object=enumeration
                        )
                        enumeration.add_conformer(conformer, auto_update=True)
                number_successful = len(
                    [
                        True
                        for enum in compound.get_enumerations()
                        if enum[0].get_molecule() is not None
                    ]
                )
                self._logger.log(
                    f"Embedding for compound {compound.get_index_string()} (name: {compound.get_name()}) completed ({number_successful} of {len(compound)} enumerations successful).",
                    _LE.INFO,
                )
            else:
                ValueError(
                    f'Value "{embed_as}" for parameter "embed_as" not supported.'
                )
