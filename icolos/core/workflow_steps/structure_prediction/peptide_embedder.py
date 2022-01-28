from pydantic import BaseModel
from icolos.core.workflow_steps.step import StepBase
from Bio import SeqIO
from Bio.PDB import PDBIO
import os
import PeptideBuilder
from PeptideBuilder import Geometry


class StepPeptideEmbedder(StepBase, BaseModel):
    def __init__(self, **data):
        super().__init__(**data)

    def execute(self):
        # use the PeptideBuilder python library to build a rough peptide structure using
        # sensible psi, phi angles etc, for subsequent simulation

        tmp_dir = self._make_tmpdir()
        self.data.generic.write_out_all_files(tmp_dir)
        # Extract the peptide sequence from the provided fasta file
        fasta_file = self.data.generic.get_argument_by_extension("fasta")

        sequences = list(SeqIO.parse(os.path.join(tmp_dir, fasta_file), format="fasta"))

        for idx, seq in enumerate(sequences):
            geom = [Geometry.geometry(aa) for aa in seq.seq]
            structure = PeptideBuilder.make_structure_from_geos(geom)

            out = PDBIO()
            out.set_structure(structure)
            # TODO: find a better naming strategy than this
            out.save(os.path.join(tmp_dir, f"sequence_{idx}.pdb"))

        self._parse_output(tmp_dir)

        self._remove_temporary(tmp_dir)
