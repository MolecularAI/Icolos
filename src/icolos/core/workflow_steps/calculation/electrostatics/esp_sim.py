from copy import deepcopy
import tempfile
from typing import List
from icolos.core.containers.compound import Conformer, Enumeration
from icolos.core.workflow_steps.step import StepBase
from pydantic import BaseModel

try:
    from espsim import GetEspSim, GetShapeSim
except ImportError:
    print(
        "WARNING - Could not import module espsim, check it is installed in your environment"
    )

from rdkit.Chem import AllChem, Mol
from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem import AllChem
from icolos.core.workflow_steps.step import _LE
from icolos.utils.general.parallelization import SubtaskContainer, Parallelizer
import os

# Based on https://github.com/hesther/espsim


class StepEspSim(StepBase, BaseModel):
    def __init__(self, **data):
        super().__init__(**data)

    def _compute_esp_sim(self, ref: Mol, trg: Enumeration, tmp_dir: str):
        """
        :param ref : Reference molecule, the known binder against which to calculate similarity
        :param trg: Icolos enumeration of the target molecule, as a smiles string.  Embedded with RDKit
        """
        # Create mol object from trg smile string
        # housekeeping for data appending later

        trg_mol = Chem.AddHs(Chem.MolFromSmiles(trg.get_smile()))
        AllChem.EmbedMolecule(ref, AllChem.ETKDG())
        AllChem.EmbedMolecule(trg_mol, AllChem.ETKDG())
        mols = [Chem.RemoveHs(x) for x in [ref, trg_mol]]
        mcs = rdFMCS.FindMCS(
            mols,
            threshold=0.8,
            completeRingsOnly=True,
            ringMatchesRingOnly=True,
        )
        self._logger.log(
            f"Computed mcs: {mcs} for enumeration {trg.get_smile()}", _LE.DEBUG
        )

        patt = Chem.MolFromSmarts(mcs.smartsString)
        refMol = mols[0]
        refMatch = refMol.GetSubstructMatch(patt)

        mv = trg_mol.GetSubstructMatch(patt)
        AllChem.AlignMol(trg_mol, refMol, atomMap=list(zip(mv, refMatch)))

        # Generate an single conformer for the reference

        # All atom alignment of the two

        charge_method = self.get_additional_setting("charge_method", default="am1-bcc")
        esp_sim = GetEspSim(trg_mol, ref, partialCharges=charge_method)
        shape_sim = GetShapeSim(trg_mol, ref)
        self._logger.log(
            f"Computed EspSim: {esp_sim}, ShapeSim: {shape_sim} for mol {trg.get_smile()}",
            _LE.DEBUG,
        )

        # now attach the mols as conformers attach the scores to the mol objects
        trg_conf = Conformer(conformer=trg_mol)
        trg_conf.get_molecule().SetProp("shape_sim", str(shape_sim))
        trg_conf.get_molecule().SetProp("esp_sim", str(esp_sim))

        trg_conf.write(os.path.join(tmp_dir, "conformer.sdf"))

    def _get_arguments(self, std_args: List) -> List:

        for flag in self.settings.arguments.flags:
            std_args.append(flag)
        for key, value in self.settings.arguments.parameters.items():
            std_args.append(f"{key}=")
            std_args.append(value)
        return std_args

    def _prepare_batch(self, batch):
        target_enums = []
        tmp_dirs = []

        for sublist in batch:
            for task in sublist:
                target_enums.append(task.data)
                tmp_dirs.append(tempfile.mkdtemp())
        return target_enums, tmp_dirs

    def _parse_output(self, trgs: List[Enumeration], tmp_dirs: List[str]) -> None:
        for tmp_dir, trg in zip(tmp_dirs, trgs):
            # grab the written sdf object
            sdf_path = os.path.join(tmp_dir, "conformer.sdf")
            mol_supplier = Chem.SDMolSupplier(sdf_path, removeHs=False)
            for mol in mol_supplier:  # should only be one conformer!
                conf = Conformer(conformer=mol)
                comp = self.get_compound_by_name(trg.get_compound_name())
                comp.find_enumeration(trg.get_enumeration_id()).add_conformer(conf)

    def _execute_espsim_parallel(self):
        # embed the reference compound

        parallelizer = Parallelizer(func=self._compute_esp_sim)
        ref_compound = Chem.AddHs(
            Chem.MolFromSmiles(self.settings.additional["ref_smiles"])
        )
        while self._subtask_container.done() is False:
            next_batch = self._get_sublists(get_first_n_lists=self._get_number_cores())

            _ = [sub.increment_tries() for element in next_batch for sub in element]
            _ = [sub.set_status_failed() for element in next_batch for sub in element]

            trgs, tmp_dirs = self._prepare_batch(next_batch)

            refs = [ref_compound for _ in range(len(next_batch))]

            parallelizer.execute_parallel(ref=refs, trg=trgs, tmp_dir=tmp_dirs)
            # hand over the embedded reference (compute once) and target compound (smiles string to be embedded)
            self._parse_output(tmp_dirs=tmp_dirs, trgs=trgs)

            for task, tmp_dir in zip(next_batch, tmp_dirs):
                for subtask in task:
                    if os.path.isfile(os.path.join(tmp_dir, "conformer.sdf")):
                        subtask.set_status_success()
                    else:
                        subtask.set_status_failed()

            self._remove_temporary(tmp_dirs)

    def execute(self):
        """
        esp-sim does molecular alignment with RDkit, then computes coulombic overlap integral + tanimoto similarity for shape measurement

        Use case takes a reference compound (known binder) and compare to REINVENT compounds

        Usage:
        * Define reference compound using settings.additional, as a smile string, to be embedded by RDkit
        * The remaining compounds are embedded using a preceeding RDkit embedding
        * attach the resulting scores to the enumeration
        """

        all_enums = []
        for compound in self.get_compounds():
            for enumeration in compound:
                all_enums.append(deepcopy(enumeration))
        if self.get_additional_setting("charge_method", default="am1-bcc") == "resp":
            # resp doesn#t play well with Icolos parallelization
            ref_compound = Chem.AddHs(
                Chem.MolFromSmiles(self.settings.additional["ref_smiles"])
            )
            tmp_dir = tempfile.mkdtemp()
            tmp_dirs = []
            for enum in all_enums:
                tmp_dir = tempfile.mkdtemp()
                tmp_dirs.append(tmp_dir)
                self._compute_esp_sim(ref_compound, enum, tmp_dir=tmp_dir)
            self._parse_output(all_enums, tmp_dirs=tmp_dirs)
            self._remove_temporary(tmp_dirs)
        else:
            self.execution.parallelization.max_length_sublists = 1
            # unroll the provided compounds,
            self._subtask_container = SubtaskContainer(max_tries=3)
            self._subtask_container.load_data(all_enums)
            self._execute_espsim_parallel()
