from pydantic import BaseModel
from rdkit import Chem

from icolos.utils.enums.program_parameters import OpenBabelEnum
from icolos.utils.enums.step_enums import StepGoldTargetPreparationEnum
from icolos.utils.execute_external.autodockvina import AutoDockVinaExecutor
from icolos.utils.execute_external.openbabel import OpenBabelExecutor
from icolos.utils.general.icolos_exceptions import StepFailed

from icolos.core.workflow_steps.step import _LE, StepBase

_STE = StepGoldTargetPreparationEnum()
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


class StepGoldTargetPreparation(StepBase, BaseModel):
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

    def specify_cavity(self):
        self._target_dict[self._TK.CAVITY_METHOD] = self._run_parameters[self._TP.CAVITY][self._TP.CAVITY_METHOD]
        if self._run_parameters[self._TP.CAVITY][self._TP.CAVITY_METHOD] == self._TP.CAVITY_METHOD_REFERENCE:
            # note, that "MoleculeReader" is able to discern many formats from the ending, including "mol2" and "pdb"
            ref_ligand_path = self._run_parameters[self._TP.CAVITY][self._TP.CAVITY_REFERENCE_PATH]
            ref_ligand = MoleculeReader(filename=ref_ligand_path)
            protein = self._settings.proteins[0]
            self._settings.binding_site = self._settings.BindingSiteFromLigand(protein, ref_ligand, distance=
            self._run_parameters[self._TP.CAVITY][self._TP.CAVITY_REFERENCE_DISTANCE])

            # add information to dictionary
            with open(ref_ligand_path, 'r') as file:
                self._target_dict[self._TK.REFERENCE_LIGAND] = [line for line in file]
            self._target_dict[self._TK.CAVITY_REFERENCE_DISTANCE] = self._run_parameters[self._TP.CAVITY][
                self._TP.CAVITY_REFERENCE_DISTANCE]
            self._target_dict[self._TK.REFERENCE_LIGAND_FILENAME] = os.path.basename(ref_ligand_path)
        elif self._run_parameters[self._TP.CAVITY][self._TP.CAVITY_METHOD] == self._TP.CAVITY_METHOD_POINT:
            raise NotImplementedError
            # origin (x,x,x)
            # distance x
        else:
            raise TargetPreparationFailed("Specified cavity determination method not defined for GOLD.")
        self._logger.log(
            f"Generated GOLD Protein.BindingSite with method {self._run_parameters[self._TP.CAVITY][self._TP.CAVITY_METHOD]}.",
            self._TL.DEBUG)

    def write_target(self, path):
        _, file_extension = os.path.splitext(path)
        if file_extension != ".pkl":
            raise TargetPreparationFailed("Receptor files must end on .pkl.")
        if self._TK.CAVITY_METHOD not in self._target_dict:
            self._logger.log("Need to have executed specify_cavity before writing out result - will attempt this now.",
                             self._TL.WARNING)
            self.specify_cavity()
        with open(path, "wb") as f:
            pickle.dump(self._target_dict, f)
        self._logger.log(f"Wrote binding site to file {path}.", self._TL.DEBUG)

    def execute(self):
        # translate input PDB file into output PDBQT file
        self._export_as_pdb2pdbqt()

        # extract and log the "box" dimensions based on the reference ligand
        self._log_extract_box()



"""import os
import pickle
import ccdc
from ccdc.docking import Docker
from ccdc.io import MoleculeReader

from dockstream.core.target_preparator import TargetPreparator

from dockstream.utils.dockstream_exceptions import TargetPreparationFailed

from dockstream.utils.enums.Gold_enums import GoldTargetPreparationEnum, GoldTargetKeywordEnum
from dockstream.containers.target_preparation_container import TargetPreparationContainer


class GoldTargetPreparator(TargetPreparator):


    def __init__(self, conf: TargetPreparationContainer, target, run_number=0):
        self._TP = GoldTargetPreparationEnum()
        self._TK = GoldTargetKeywordEnum()
        self._target_dict = {self._TK.VERSION: self._TK.CURRENT_VERSION}

        # invoke base class's constructor first
        super().__init__(conf=conf, run_number=run_number)

        # check, whether the backend run specified is an "GOLD" one
        if self._run_parameters[self._TP.RUNS_BACKEND] != self._TP.RUNS_BACKEND_GOLD:
            raise TargetPreparationFailed("Tried to make an GOLD preparation with different backend specification.")

        if isinstance(target, str):
            if os.path.isfile(target):
                _, file_extension = os.path.splitext(target)
                if file_extension == ".pdb":
                    self._target = Docker()
                    self._settings = self._target.settings
                    self._settings.add_protein_file(file_name=target)

                    # add information to dictionary
                    with open(target, 'r') as file:
                        self._target_dict[self._TK.TARGET_PDB] = [line for line in file]
                    self._target_dict[self._TK.TARGET_PDB_FILENAME] = os.path.basename(target)
                else:
                    raise TargetPreparationFailed("Specified input file must be in PDB format for GOLD.")
            else:
                raise TargetPreparationFailed("Input target file does not exist.")
        elif isinstance(target, ccdc.docking.Docker):
            raise NotImplementedError
        else:
            raise TargetPreparationFailed("Constructor only accepts and Protein.BindingSite object or a file path.")
        self._logger.log("Added target to GOLD settings.", self._TL.DEBUG)
"""