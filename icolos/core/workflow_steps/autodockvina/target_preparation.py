from pydantic import BaseModel
from rdkit import Chem

from icolos.utils.enums.program_parameters import OpenBabelEnum
from icolos.utils.enums.step_enums import StepAutoDockVinaTargetPreparationEnum
from icolos.utils.execute_external.autodockvina import AutoDockVinaExecutor
from icolos.utils.execute_external.openbabel import OpenBabelExecutor
from icolos.utils.general.icolos_exceptions import StepFailed

from icolos.core.workflow_steps.step import _LE, StepBase

_STE = StepAutoDockVinaTargetPreparationEnum()
_OBE = OpenBabelEnum()


class ADVExtractBoxTP(BaseModel):
    reference_ligand_path: str = None
    reference_ligand_format: str = _STE.EXTRACT_BOX_REFERENCE_LIGAND_FORMAT_PDB


class ADVAdditionalTP(BaseModel):
    pH: float = (
        7.4  # set target pH value that determines the protein's side-chain states
    )
    input_receptor_pdb: str = None
    output_receptor_pdbqt: str = None
    extract_box: ADVExtractBoxTP = ADVExtractBoxTP()


class StepAutoDockVinaTargetPreparation(StepBase, BaseModel):
    _openbabel_executor: OpenBabelExecutor = None
    adv_additional: ADVAdditionalTP = None

    def __init__(self, **data):
        super().__init__(**data)

        # initialize the executor and test availability
        self._initialize_backend(executor=AutoDockVinaExecutor)
        self._check_backend_availability()

        # initialize the executor for all "OpenBabel"
        self._openbabel_executor = OpenBabelExecutor()
        if not self._openbabel_executor.is_available():
            raise StepFailed(
                "AutoDock Vina requires OpenBabel execution, initialization failed."
            )

        # set ADV specific settings and ensure that each molecule gets its own sublist
        self.adv_additional = ADVAdditionalTP(**self.settings.additional)

    def _export_as_pdb2pdbqt(self):
        # Note: In contrast to the ligand preparation, we will not use a tree-based flexibility treatment here - thus,
        #       the option "-xr" is used. Partial charges of the receptor are not used in AutoDock Vina.
        arguments = [
            " ".join(
                [_OBE.OBABEL_INPUTFORMAT_PDB, self.adv_additional.input_receptor_pdb]
            ),
            _OBE.OBABEL_OUTPUT_FORMAT_PDBQT,
            " ".join([_OBE.OBABEL_O, self.adv_additional.output_receptor_pdbqt]),
            "".join([_OBE.OBABEL_X, _OBE.OBABEL_X_R]),
            _OBE.OBABEL_P,
            str(self.adv_additional.pH),
            _OBE.OBABEL_PARTIALCHARGE,
            _OBE.OBABEL_PARTIALCHARGE_GASTEIGER,
        ]
        self._openbabel_executor.execute(
            command=_OBE.OBABEL, arguments=arguments, check=True
        )
        self._logger.log(
            f"Exported target as PDBQT file {self.adv_additional.output_receptor_pdbqt}.",
            _LE.INFO,
        )

    def _log_extract_box(self):
        x_coords, y_coords, z_coords = self._extract_box()
        if x_coords is not None:

            def dig(value):
                return round(value, ndigits=2)

            self._logger.log(
                f"Calculating lingad dimensions for AutoDock Vina docking protocol.",
                _LE.INFO,
            )
            self._logger.log(
                f"Ligand ({self.adv_additional.extract_box.reference_ligand_path}):",
                _LE.INFO,
            )
            self._logger_blank.log(
                f"X coordinates: min={dig(min(x_coords))}, max={dig(max(x_coords))}, mean={dig(sum(x_coords) / len(x_coords))}",
                _LE.INFO,
            )
            self._logger_blank.log(
                f"Y coordinates: min={dig(min(y_coords))}, max={dig(max(y_coords))}, mean={dig(sum(y_coords) / len(y_coords))}",
                _LE.INFO,
            )
            self._logger_blank.log(
                f"Z coordinates: min={dig(min(z_coords))}, max={dig(max(z_coords))}, mean={dig(sum(z_coords) / len(z_coords))}",
                _LE.INFO,
            )

    def _extract_box(self):
        # extracts box suggestions from a reference ligand, which can be added to a AutoDock Vina run
        # load the reference file (PDB or SDF)
        ref_format = self.adv_additional.extract_box.reference_ligand_format.upper()
        if ref_format == _STE.EXTRACT_BOX_REFERENCE_LIGAND_FORMAT_PDB:
            ref_mol = Chem.MolFromPDBFile(
                self.adv_additional.extract_box.reference_ligand_path, sanitize=True
            )
        elif ref_format == _STE.EXTRACT_BOX_REFERENCE_LIGAND_FORMAT_SDF:
            mol_supplier = Chem.SDMolSupplier(
                self.adv_additional.extract_box.reference_ligand_path
            )
            for mol in mol_supplier:
                if mol is None:
                    raise StepFailed(
                        f"Could not load molecule from {self.adv_additional.extract_box.reference_ligand_path} - abort."
                    )
                ref_mol = mol
                break
        else:
            raise StepFailed(
                f"Reference ligand format {ref_format} not supported, use PDB or SDF instead - abort."
            )

        # extract coordinates
        x_coords = [atom[0] for atom in ref_mol.GetConformer(0).GetPositions()]
        y_coords = [atom[1] for atom in ref_mol.GetConformer(0).GetPositions()]
        z_coords = [atom[2] for atom in ref_mol.GetConformer(0).GetPositions()]
        return x_coords, y_coords, z_coords

    def execute(self):
        # translate input PDB file into output PDBQT file
        self._export_as_pdb2pdbqt()

        # extract and log the "box" dimensions based on the reference ligand
        self._log_extract_box()
