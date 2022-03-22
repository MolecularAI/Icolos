from copy import deepcopy
from icolos.core.containers.generic import GenericContainer, GenericData
import json
import pandas as pd
from rdkit import Chem
from icolos.core.containers.gmx_state import GromacsState

from icolos.loggers.base_logger import BaseLogger
from icolos.utils.enums.input_enums import InputEnum
from icolos.utils.enums.logging_enums import LoggingConfigEnum
from icolos.utils.enums.write_out_enums import WriteOutEnum
from icolos.utils.general.icolos_exceptions import StepFailed
from icolos.utils.smiles import to_smiles
from icolos.utils.general.files_paths import infer_input_type

from typing import List, Any
from pydantic import BaseModel

from icolos.core.step_utils.input_merger import InputMerger, StepMerge
from icolos.core.containers.compound import Enumeration, Compound, Conformer
from icolos.utils.enums.step_enums import StepBaseEnum, StepGromacsEnum
import os
from tempfile import mkdtemp
import requests

_SBE = StepBaseEnum
_LE = LoggingConfigEnum()
_WE = WriteOutEnum()
_SGE = StepGromacsEnum()
_IE = InputEnum()


class StringPath(str):
    def __new__(cls, content):
        return super().__new__(cls, content)


class StringFile(str):
    def __new__(cls, content):
        return super().__new__(cls, content)


class StepData(BaseModel):
    class Config:
        arbitrary_types_allowed = True

    compounds: List[Compound] = []
    generic: GenericContainer = GenericContainer()
    gmx_state: GromacsState = GromacsState()


class StepCSVInputColumnParameters(BaseModel):
    smiles: str
    names: str = None


class StepInputEnforceIDs(BaseModel):
    compound_ids: List = None
    enumeration_ids: List = None


class StepInputSource(BaseModel):
    source: str
    source_type: str = None
    source_field: str = _IE.SOURCE_FIELD_COMPOUNDS
    target_field: str = _IE.SOURCE_FIELD_COMPOUNDS
    extension: str = None
    format: str = None
    delimiter: str = ","
    columns: StepCSVInputColumnParameters = None
    enforce_ids: StepInputEnforceIDs = None


class StepInputParameters(BaseModel):

    compounds: List[StepInputSource] = []
    generic: List[StepInputSource] = []
    gmx_state: StepInputSource = {}
    perturbation_map: List[StepInputSource] = None
    merge: StepMerge = StepMerge()
    work_dir: str = None


class InputPreparator(BaseModel):
    workflow: Any = None
    logger: BaseLogger = None

    class Config:
        underscore_attrs_are_private = True
        arbitrary_types_allowed = True

    def __init__(self, **data):
        super().__init__(**data)

    def generate_input(self, step_input: StepInputParameters, step_type):
        compounds = self._generate_compound_input(step_input)
        generic = (
            self._generate_generic_input(step_input, step_type)
            if step_input.generic
            else GenericContainer()
        )
        gmx_state = (
            self._generate_gmx_state_input(step_input)
            if step_input.gmx_state
            else GromacsState()
        )
        # Instruct the step to run in a specific workdir, e.g. from a previously failed job or to execute a few related steps in the same dir
        if step_input.work_dir is not None:
            if os.path.isdir(step_input.work_dir):
                work_dir = step_input.work_dir
                self.logger.log(
                    f"Found specified work dir at {step_input.work_dir}", _LE.DEBUG
                )
                # now check whether this needs attaching to the workflow for the rest of the steps
                if self.workflow.header.global_settings.single_directory:
                    self.workflow.workflow_data.work_dir = work_dir
                    self.logger.log(
                        f"Setting workdir  at {step_input.work_dir} to the workflow's workdir",
                        _LE.DEBUG,
                    )
            else:
                # last resort, if a previous step_id has been passed, get the work_dir from there

                work_dir = self.workflow.find_step_by_step_id(
                    step_input.work_dir
                ).work_dir
        elif (
            self.workflow is not None
            and self.workflow.header.global_settings.single_directory
        ):
            # Entire workflow running in a single dir (pmx), either generate one for the first
            # step or use the already generated dir
            work_dir = self._get_workflow_workdir()
        else:
            work_dir = None
        return (
            StepData(compounds=compounds, generic=generic, gmx_state=gmx_state),
            work_dir,
        )

    def _get_workflow_workdir(self):
        # check whether the workflow already has one attached, otherwise create one
        if self.workflow.workflow_data.work_dir is not None and os.path.isdir(
            self.workflow.workflow_data.work_dir
        ):
            return self.workflow.workflow_data.work_dir
        else:
            tmp_dir = mkdtemp()
            self.workflow.workflow_data.work_dir = tmp_dir
            self.logger.log(f"Set workflow's tmpdir to {tmp_dir}", _LE.DEBUG)
            return tmp_dir

    def _generate_compound_input(self, step_input: StepInputParameters) -> List:
        compounds = []
        for inp in step_input.compounds:
            if inp.target_field == _IE.TARGET_FIELD_COMPOUNDS:
                buffer = []
                if inp.source_type == _SBE.INPUT_SOURCE_TYPE_FILE:
                    buffer.append(self._read_compound_input_from_file(inp))
                elif inp.source_type == _SBE.INPUT_SOURCE_TYPE_STEP:
                    buffer.append(self._read_compound_input_from_step(inp))
                elif inp.source_type == _SBE.INPUT_SOURCE_TYPE_STRING:
                    buffer.append(self._read_input_from_string(inp))
                else:
                    raise ValueError(
                        f"Source type {inp.source_type} for compound input unsupported - abort."
                    )
                if inp.target_field == _IE.SOURCE_FIELD_COMPOUNDS:
                    # note: no unrolling here!
                    compounds = compounds + buffer

            elif inp.target_field == _IE.TARGET_FIELD_CONFORMERS:
                if inp.source_type == _SBE.INPUT_SOURCE_TYPE_FILE:
                    compounds = compounds + self._read_conformers_input_from_file(inp)
        if len(compounds) > 0:
            compounds = self._apply_compound_merger(step_input, compounds)
        return compounds

    def _generate_gmx_state_input(
        self, step_input: StepInputParameters
    ) -> GromacsState:
        input_step = self.workflow.find_step_by_step_id(step_input.gmx_state.source)
        return deepcopy(input_step.data.gmx_state)

    def _generate_generic_input(
        self, step_input: StepInputParameters, step_type
    ) -> GenericContainer:
        generic = GenericContainer()
        for inp in step_input.generic:
            files = self._read_data_to_generic(inp)
            generic.add_files(files)
        return generic

    def _read_data_to_generic(self, inp: StepInputSource):
        ext = inp.extension
        if inp.source_type == _SBE.INPUT_SOURCE_TYPE_FILE or os.path.isfile(inp.source):
            assert os.path.isfile(inp.source)
            try:
                with open(inp.source, "r") as f:
                    data = f.read()
            except UnicodeDecodeError:
                with open(inp.source, "rb") as f:
                    data = f.read()
            file = GenericData(inp.source.split("/")[-1], data)
            return [file]
        elif inp.source_type == _SBE.INPUT_SOURCE_TYPE_URL or inp.source.startswith(
            "http"
        ):
            data = self._get_pdb_file_from_api(inp.source)
            file_name = inp.source.split("/")[-1].split(".")[0] + "." + inp.extension
            file = GenericData(file_name=file_name, file_data=data)
            return [file]
        elif inp.source_type == _SBE.INPUT_SOURCE_TYPE_DIR or os.path.isdir(inp.source):
            assert os.path.isdir(inp.source)
            file = GenericData(
                file_data=inp.source,
                file_name=inp.source.split("/")[-1],
                extension=inp.extension,
            )
            return [file]
        else:
            # fall back on step source type
            input_step = self.workflow.find_step_by_step_id(inp.source)
            files = input_step.data.generic.get_files_by_extension(ext)

            # special case for itp and ndx files, these are included in the topol file so are never arguments
            if ext in ["itp", "ndx"]:
                return files

            if len(files) == 1:
                file = files[0]
                return [file]
            # else use the argument method
            else:
                # this introduces a manual check on which file the user wants if there are multiple
                file = input_step.data.generic.get_argument_by_extension(
                    ext, rtn_file_object=True
                )
                return [file]

    def _get_pdb_file_from_api(self, pdb_url: str):
        response = self._get_request(pdb_url)
        if response is None or not response.ok:
            return None
        return response.text

    def _get_request(self, url, max_tries=5):
        trials = 0
        while trials < max_tries:
            response = requests.get(url)
            if response.status_code == 200:
                return response

    def _apply_compound_merger(
        self, step_input: StepInputParameters, compounds: List[Compound]
    ) -> List[Compound]:
        merger = InputMerger(step_input.merge)
        compounds = merger.unroll_compounds(compounds)
        if not any(
            [
                True
                for compound in step_input.compounds
                if compound.enforce_ids is not None
            ]
        ):
            compounds = merger.merge(compounds=compounds)

        if len(compounds) == 0 and self.logger is not None:
            self.logger.log(
                "Input list of compounds is empty, this is likely an error.",
                _LE.WARNING,
            )
        return compounds

    def _read_compound_input_from_step(self, inp: StepInputSource):
        input_step = self.workflow.find_step_by_step_id(inp.source)
        return input_step.clone_compounds()

    def _read_conformers_input_from_file(self, inp: StepInputSource):
        # set up path to input file and extract the input format
        input_format = inp.format
        if input_format is None and self.logger is not None:
            self.logger.log(
                "No input format specified, will try to infer type (not recommended).",
                _LE.WARNING,
            )
            input_format = infer_input_type(inp.source)
        input_format = input_format.upper()

        # call the respective loading function
        if input_format == _SBE.FORMAT_SDF:
            compound = Compound(compound_number=0)
            enumeration = Enumeration()
            for mol_id, mol in enumerate(
                Chem.SDMolSupplier(inp.source, removeHs=False)
            ):
                conformer = Conformer(conformer=mol, enumeration_object=enumeration)
                enumeration.add_conformer(conformer=conformer, auto_update=True)
            compound.add_enumeration(enumeration, auto_update=True)
            return [compound]
        else:
            raise ValueError(
                f"At the moment, input format {input_format} is not supported."
            )

    def _read_compound_input_from_file(self, inp: StepInputSource):
        # set up path to input file and extract the input format
        input_format = inp.format
        if input_format is None and self.logger is not None:
            self.logger.log(
                "No input format specified, will try to infer type (not recommended).",
                _LE.WARNING,
            )
            input_format = infer_input_type(inp.source)
        input_format = input_format.upper()

        # call the respective loading function
        if input_format == _SBE.FORMAT_SDF:
            result = self._read_in_SDF_file(inp)
        elif input_format == _SBE.FORMAT_CSV:
            result = self._read_in_CSV_file(inp)
        elif input_format == _SBE.FORMAT_SMI:
            result = self._read_in_SMI_file(inp)
        elif input_format == _SBE.FORMAT_JSON:
            result = self._read_in_JSON_file(inp)
        else:
            raise ValueError(
                f"At the moment, input format {input_format} is not supported."
            )

        # apply ID enforcement, if specified
        return self._enforce_ids(result, inp)

    def _read_input_from_string(self, inp: StepInputSource) -> List[Compound]:
        # the strings must be separated by a semi-colon (';'); they may have names in front separated by a colon (':')
        elements = inp.source.split(";")
        list_compounds = []
        for line_id, line in enumerate(elements):
            # it could be, that names are part of the elements; otherwise use the number
            # remove trailing or preceding white spaces
            parts = [x.strip() for x in line.split(":")]
            if len(parts) == 2:
                compound = Compound(name=parts[0], compound_number=line_id)
                smile = parts[1]
            else:
                compound = Compound(name=str(line_id), compound_number=line_id)
                smile = parts[0]
            enumeration = Enumeration(smile=smile, original_smile=smile)
            compound.add_enumeration(enumeration, auto_update=True)
            list_compounds.append(compound)

        # apply ID enforcement, if specified
        return self._enforce_ids(list_compounds, inp)

    def _read_in_SDF_file(self, inp: StepInputSource) -> List[Compound]:
        def _get_existing_enumeration(comp_id, enum_id):
            comp = _get_existing_compound(comp_id)
            for enum in comp.get_enumerations():
                if enum.get_enumeration_id() == int(enum_id):
                    return enum
            raise ValueError

        def _get_existing_compound(idx):
            for comp in list_compounds:
                if int(idx) == comp.get_compound_number():
                    return comp
            raise ValueError

        list_compounds = []
        compound_number = 0
        icolos_naming = True
        # Parses compounds following the Icolos naming convention of Compound:Enumeration:Conformer to reconstruct the compound object
        for mol in Chem.SDMolSupplier(inp.source, removeHs=False):
            new_compound = False
            new_enumeration = False
            mol_name = mol.GetProp(_WE.RDKIT_NAME)
            # assuming the mol name follows Icolos conventions
            try:
                id_parts = mol_name.split(":")
                comp_id = id_parts[0]
                enum_id = id_parts[1]

            except:
                #
                icolos_naming = False
                comp_id = mol_name
                enum_id = 0

            if icolos_naming:
                # reconstruct compound objects
                try:
                    # try to find an existing compound with the correct name
                    compound = _get_existing_compound(idx=comp_id)
                except ValueError:
                    # the compound does not yet exist, create the object
                    new_compound = True
                    try:
                        # if we have standard icolos compound naming
                        comp_num = int(comp_id)
                    except ValueError:
                        # some other naming scheme
                        comp_num = compound_number
                    compound = Compound(name=comp_id, compound_number=comp_num)
                try:
                    # check whether the enumeration exists
                    enumeration = _get_existing_enumeration(comp_id, enum_id)
                except ValueError:
                    new_enumeration = True
                    enumeration = Enumeration(
                        smile=to_smiles(mol),
                        molecule=mol,
                        original_smile=to_smiles(mol),
                    )

                if len(id_parts) == 3:
                    # i.e. 0:0:0, we have a conformer
                    conf = Conformer(
                        conformer=mol,
                        enumeration_object=enumeration,
                        conformer_id=int(id_parts[2]),
                    )
                    enumeration.add_conformer(conf, auto_update=True)
                if new_enumeration:
                    compound.add_enumeration(enumeration, auto_update=True)
                if new_compound:
                    list_compounds.append(compound)

            else:
                # if non-standard naming conventions, simply load each mol into a new compound object, with single enum/conf
                compound = Compound(name=comp_id, compound_number=compound_number)
                enum = Enumeration(
                    smile=to_smiles(mol),
                    molecule=mol,
                    original_smile=to_smiles(mol),
                    enumeration_id=0,
                )
                enum.add_conformer(
                    Conformer(conformer=mol, enumeration_object=mol, conformer_id=0),
                    auto_update=True,
                )
                compound.add_enumeration(enumeration=enum)
                list_compounds.append(compound)

            compound_number += 1
        return list_compounds

    def _read_in_SMI_file(self, inp: StepInputSource) -> List[Compound]:
        list_compounds = []
        with open(inp.source, "r") as f:
            # while the SMI file definition requires a name (separated by blanks) for each line
            # as well, assume that this might not be present
            lines = [line.rstrip() for line in f.readlines()]
        for line_id, line in enumerate(lines):
            if line == "":
                continue

            parts = line.split()
            if len(parts) == 2:
                compound = Compound(name=parts[1], compound_number=line_id)
            else:
                compound = Compound(name=str(line_id), compound_number=line_id)
            enumeration = Enumeration(smile=parts[0], original_smile=parts[0])
            compound.add_enumeration(enumeration, auto_update=True)
            list_compounds.append(compound)
        return list_compounds

    def _read_in_JSON_file(self, inp: StepInputSource) -> List[Compound]:
        list_compounds = []

        # load input
        with open(inp.source, "r") as f:
            inp_json = f.read().replace("\r", "").replace("\n", "")
            inp_dict = json.loads(inp_json)

        comp_id = 0
        for name, smile in zip(inp_dict[_IE.JSON_NAMES], inp_dict[_IE.JSON_SMILES]):
            compound = Compound(name=name, compound_number=comp_id)
            enumeration = Enumeration(smile=smile, original_smile=smile)
            compound.add_enumeration(enumeration, auto_update=True)
            list_compounds.append(compound)
            comp_id += 1

        return list_compounds

    def _read_in_CSV_file(self, inp: StepInputSource) -> List[Compound]:
        list_compounds = []
        delimiter = inp.delimiter
        data = pd.read_csv(inp.source, delimiter=delimiter)

        smiles_column = inp.columns.smiles
        if smiles_column not in list(data.columns):
            raise StepFailed(
                f"Column name for the smiles either not set or not found in input CSV."
            )

        # deal with names (if specified)
        names_column = inp.columns.names
        if names_column is None:
            names_compounds = None
        else:
            if names_column not in list(data.columns):
                raise StepFailed(
                    f"Specified column name ({names_column}) for the names either not found in input CSV."
                )
            else:
                names_compounds = [
                    str(name).strip() for name in data[names_column].tolist()
                ]

        # build the compounds
        smiles = [str(line).strip() for line in data[smiles_column].tolist()]
        for number in range(len(smiles)):
            if names_compounds is not None:
                compound = Compound(
                    name=names_compounds[number], compound_number=number
                )
            else:
                compound = Compound(name=str(number), compound_number=number)
            enumeration = Enumeration(
                smile=smiles[number], original_smile=smiles[number]
            )
            compound.add_enumeration(enumeration, auto_update=True)
            list_compounds.append(compound)
        return list_compounds

    def _enforce_ids(
        self, compounds: List[Compound], inp: StepInputSource
    ) -> List[Compound]:
        if inp.enforce_ids is not None:
            if inp.enforce_ids.compound_ids is not None:
                for comp_idx, comp in enumerate(compounds):
                    comp.set_compound_number(
                        int(inp.enforce_ids.compound_ids[comp_idx])
                    )

            # set enumeration ids
            enum_id_idx = 0
            if inp.enforce_ids.enumeration_ids is not None:
                for comp in compounds:
                    for enum in comp.get_enumerations():
                        enum.set_enumeration_id(
                            int(inp.enforce_ids.enumeration_ids[enum_id_idx])
                        )
                        enum_id_idx += 1
            if self.logger is not None:
                self.logger.log(
                    "Enforced IDs for compounds and enumerations specified (merging disabled).",
                    _LE.DEBUG,
                )
        return compounds
