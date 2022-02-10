import os
from typing import Optional, Iterable, Union
from pydantic import BaseModel
from icolos.core.workflow_steps.step import StepBase
from icolos.utils.enums.program_parameters import SchrodingerExecutablesEnum
from icolos.utils.execute_external.license_token_guard import (
    TokenGuardParameters,
    SchrodingerLicenseTokenGuard,
)
from icolos.utils.general.files_paths import attach_root_path
from icolos.utils.execute_external.schrodinger import SchrodingerExecutor
from icolos.utils.enums.step_enums import StepDesmondEnum
from icolos.core.workflow_steps.step import _LE
import re
from shutil import copy
from typing import Dict


_EE = SchrodingerExecutablesEnum()
_SDE = StepDesmondEnum()


class StepSchrodingerBase(StepBase, BaseModel):

    token_guard: Optional[TokenGuardParameters] = None

    def __init__(self, **data):
        super().__init__(**data)

    def _apply_token_guard(self):
        if self.token_guard is not None:
            token_guard = SchrodingerLicenseTokenGuard(token_guard=self.token_guard)
            token_guard.guard()

    # TODO: Deprecated - use self.converter
    def _translate_SDF_to_MAE(
        self, sdf_path: str, mae_path: str, executor: SchrodingerExecutor
    ):
        """As "Glide" is only able to read MAE (Maestro) files, write the ligands out in that format."""

        # call "sdconvert" from Schrodinger's software
        arguments = [
            "".join([_EE.SDCONVERT_I, _EE.SDCONVERT_FORMAT_SD]),
            sdf_path,
            "".join([_EE.SDCONVERT_O, _EE.SDCONVERT_FORMAT_MAE]),
            mae_path,
        ]
        execution_result = executor.execute(
            command=_EE.SDCONVERT, arguments=arguments, check=True
        )

    def _translate_MAE_to_SDF(
        self, mae_path: str, sdf_path: str, executor: SchrodingerExecutor
    ):
        """In cases where the write-out mode for Glide is not producing SDF files."""

        # call "sdconvert" from Schrodinger's software
        arguments = [
            "".join([_EE.SDCONVERT_I, _EE.SDCONVERT_FORMAT_MAE]),
            mae_path,
            "".join([_EE.SDCONVERT_O, _EE.SDCONVERT_FORMAT_SD]),
            sdf_path,
        ]
        execution_result = executor.execute(
            command=_EE.SDCONVERT, arguments=arguments, check=True
        )

    def _translate_PDB_to_MAE(
        self, pdb_path: str, mae_path: str, executor: SchrodingerExecutor
    ):
        """In cases where the write-out mode for Glide is not producing SDF files."""

        # call "sdconvert" from Schrodinger's software
        arguments = [
            "".join([_EE.SDCONVERT_I, _EE.STRUCTCAT_FORMAT_PDB]),
            pdb_path,
            "".join([_EE.SDCONVERT_O, _EE.SDCONVERT_FORMAT_MAE]),
            mae_path,
        ]
        execution_result = executor.execute(
            command=_EE.STRUCTCONVERT, arguments=arguments, check=True
        )

    def _replace_config_value(self, key, value, config):
        value = str(value)
        pattern = rf"({key} =).*"
        pattern = re.compile(pattern)
        config = re.sub(pattern, rf"\1 {value}", config)
        return config

    def _get_template(self, file_name):
        file = [
            file
            for file in os.listdir(attach_root_path("icolos/config/desmond"))
            if file_name in file
        ]
        assert len(file) == 1
        return file[0]

    def _write_config(self, tmp_dir, dict_: Dict, file_name):
        # see if a config file was specified, assume no further changes:
        if _SDE.CONFIG in dict_.keys() and dict_[_SDE.CONFIG] is not None:
            copy(dict_[_SDE.CONFIG], tmp_dir)
        else:
            template = self._get_template(file_name)
            with open(attach_root_path(f"icolos/config/desmond/{template}"), "r") as f:
                config = f.read()
            for k, v in dict_.items():
                config = self._replace_config_value(k, v, config)

            self._logger.log(f"Compiled  file {file_name}...", _LE.DEBUG)
            for line in config.split("\n"):
                self._logger_blank.log(line, _LE.DEBUG)
            with open(os.path.join(tmp_dir, file_name), "w") as f:
                f.write(config)

    def _parse_arguments(self, defaults):
        args = []

        for flag in self.settings.arguments.flags:
            args.append(flag)
        if "-WAIT" not in args:
            args.append("-WAIT")
        for k, v in self.settings.arguments.parameters.items():
            args.append(k)
            args.append(v)
        for k, v in defaults.items():
            if k not in args:
                args.append(k)
                args.append(v)
        return args

    @staticmethod
    def _parse_maestro_in_file(
        lines: Iterable[str],
    ) -> Dict[str, Union[str, Dict[str, str]]]:
        """Parses Maestro input, and returns keywords dict for it."""

        separator3 = "   "
        indent4 = "    "
        block_starters = {
            "[CONSTRAINT_GROUP",
            "[FEATURE",
        }

        # All Glide keywords. Get all keywords with:
        #   $ module load schrodinger
        #   $ glide -docking-keywords | cut -d' ' -f1 | sed 's/.*/"&"/' | paste -sd , -
        # List keywords, get first word, wrap in quotes, join lines.
        # See:
        #   - https://stackoverflow.com/a/19145499
        #   - https://unix.stackexchange.com/a/251362
        allowed_keywords = {
            "AMIDE_MODE",
            "AMIDE_TRANS_ALL",
            "AMIDE_TRANSTOL",
            "ASL_RES_INTERACTION",
            "CALC_INPUT_RMS",
            "CANONICALIZE",
            "COMPRESS_POSES",
            "CORE_ATOMS",
            "CORE_DEFINITION",
            "CORE_FILTER",
            "CORE_POS_MAX_RMSD",
            "CORE_RESTRAIN",
            "CORE_RESTRAIN_V",
            "CORE_SMARTS",
            "CORE_SNAP",
            "CORECONS_FALLBACK",
            "CSV_PROPS_FILE",
            "CV_CUTOFF",
            "DIELMOD",
            "DOCKING_METHOD",
            "DOINTRA",
            "DOINTRA_SCALE",
            "DSCORE_CUTOFF",
            "EPIK_PENALTIES",
            "EXPANDED_SAMPLING",
            "FITDEN",
            "FORCEFIELD",
            "FORCEPLANAR",
            "GLIDE_CONFGEN_BADDIST2",
            "GLIDE_CONFGEN_EFCUT",
            "GLIDE_CONS_FEAT_FILE",
            "GLIDE_CONS_FINALONLY",
            "GLIDE_CONS_RMETCOORD",
            "GLIDE_CONS_RNOEMAX",
            "GLIDE_CONS_RNOEMIN",
            "GLIDE_CONS_RPOS",
            "GLIDE_CONS_XMETCOORD",
            "GLIDE_CONS_XNOE",
            "GLIDE_CONS_XPOS",
            "GLIDE_CONS_YMETCOORD",
            "GLIDE_CONS_YNOE",
            "GLIDE_CONS_YPOS",
            "GLIDE_CONS_ZMETCOORD",
            "GLIDE_CONS_ZNOE",
            "GLIDE_CONS_ZPOS",
            "GLIDE_DIELCO",
            "GLIDE_ELEMENTS",
            "GLIDE_EXVOL_PENAL_NUM",
            "GLIDE_EXVOL_PENAL_STRENGTH",
            "GLIDE_NTOTALCONS",
            "GLIDE_NUMEXVOL",
            "GLIDE_NUMMETCOORDCONS",
            "GLIDE_NUMMETCOORDSITES",
            "GLIDE_NUMNOECONS",
            "GLIDE_NUMPOSITCONS",
            "GLIDE_NUMUSEXVOL",
            "GLIDE_OUTPUT_USEHTOR",
            "GLIDE_REFLIG_FORMAT",
            "GLIDE_REXVOL",
            "GLIDE_REXVOLIN",
            "GLIDE_TORCONS_ALLBONDS",
            "GLIDE_TORCONS_IATOMS",
            "GLIDE_TORCONS_JATOMS",
            "GLIDE_TORCONS_KATOMS",
            "GLIDE_TORCONS_LATOMS",
            "GLIDE_TORCONS_PATTERN_INDEX",
            "GLIDE_TORCONS_PATTERNS",
            "GLIDE_TORCONS_SETVAL",
            "GLIDE_TORCONS_VALUES",
            "GLIDE_TORCONSFILE",
            "GLIDE_XEXVOL",
            "GLIDE_XP_NMAXCORE",
            "GLIDE_XP_RMSCUT",
            "GLIDE_YEXVOL",
            "GLIDE_ZEXVOL",
            "GLIDECONS",
            "GLIDECONSFEATATOMS",
            "GLIDECONSFEATHASINCLUDE",
            "GLIDECONSFEATINCLUDE",
            "GLIDECONSFEATINDEX",
            "GLIDECONSFEATPATTERNS",
            "GLIDECONSGROUPNREQUIRED",
            "GLIDECONSNAMES",
            "GLIDECONSUSEMET",
            "GLIDESCORUSEMET",
            "GLIDEUSEALLEXVOL",
            "GLIDEUSECONSFEAT",
            "GLIDEUSECONSFEATINDEX",
            "GLIDEUSECONSGROUPINDEX",
            "GLIDEUSECONSLABELS",
            "GLIDEUSEXVOL",
            "GLIDEUSEXVOLNAMES",
            "GLIDEXVOLNAMES",
            "GRIDFILE",
            "GSCORE",
            "GSCORE_CUTOFF",
            "HAVEGLIDECONSFEAT",
            "HBOND_ACCEP_HALO",
            "HBOND_CUTOFF",
            "HBOND_DONOR_AROMH",
            "HBOND_DONOR_AROMH_CHARGE",
            "HBOND_DONOR_HALO",
            "INCLUDE_INPUT_CONF",
            "INCLUDE_INPUT_RINGS",
            "JOBNAME",
            "KEEP_SUBJOB_POSES",
            "KEEPRAW",
            "KEEPSKIPPED",
            "LIG_CCUT",
            "LIG_MAECHARGES",
            "LIG_VSCALE",
            "LIGAND_END",
            "LIGAND_START",
            "LIGANDFILE",
            "LIGANDFILES",
            "LIGFORMAT",
            "LIGPREP",
            "LIGPREP_ARGS",
            "MACROCYCLE",
            "MACROCYCLE_OPTIONS",
            "MAX_ITERATIONS",
            "MAXATOMS",
            "MAXKEEP",
            "MAXREF",
            "MAXROTBONDS",
            "METAL_CUTOFF",
            "NENHANCED_SAMPLING",
            "NMAXRMSSYM",
            "NOSORT",
            "NREPORT",
            "NREQUIRED_CONS",
            "OUTPUTDIR",
            "PAIRDISTANCES",
            "PEPTIDE",
            "PHASE_DB",
            "PHASE_NCONFS",
            "PHASE_SUBSET",
            "POSE_DISPLACEMENT",
            "POSE_HTORSION",
            "POSE_OUTTYPE",
            "POSE_RMSD",
            "POSES_PER_LIG",
            "POSTDOCK",
            "POSTDOCK_ITMAX",
            "POSTDOCK_NPOSE",
            "POSTDOCK_SCITMAX",
            "POSTDOCK_XP_DELE",
            "POSTDOCKCG",
            "POSTDOCKLIGMIN",
            "POSTDOCKSTRAIN",
            "PRECISION",
            "PREMIN",
            "PREMINCG",
            "PREMINELEC",
            "PREMINITMAX",
            "RADIUS_RES_INTERACTION",
            "REF_LIGAND_FILE",
            "REFINDEX",
            "REPORT_CPU_TIME",
            "REWARD_INTRA_HBONDS",
            "RINGCONFCUT",
            "RINGONFLY",
            "SAMPLE_N_INVERSIONS",
            "SAMPLE_RINGS",
            "SCORE_INPUT_POSE",
            "SCORE_MINIMIZED_INPUT_POSE",
            "SCORING_CUTOFF",
            "SHAPE_ATOMS",
            "SHAPE_RESTRAIN",
            "SHAPE_TYPING",
            "SKIP_EPIK_METAL_ONLY",
            "STRAIN_GSFACTOR",
            "STRAIN_GSTHRESH",
            "STRAINELEC",
            "SUBSTRATE_PENAL_FILE",
            "USE_CONS",
            "USE_REF_LIGAND",
            "USECOMPMAE",
            "WRITE_CSV",
            "WRITE_RES_INTERACTION",
            "WRITE_TIMINGS_CSV",
            "WRITE_XP_DESC",
            "WRITEREPT",
        }

        result = {}
        current_block = None
        for linenum, line in enumerate(lines):
            if any(line.startswith(starter) for starter in block_starters):
                # Block start.
                current_block = line.strip()
                result[current_block] = {}
            elif line.strip() == "":
                # Empty line: close current block if any is open, and skip the line.
                current_block = None
            elif line.startswith(indent4):
                # Indented line inside the block.
                if current_block is None:
                    raise ValueError(
                        f"Unexpected indent outside of block for line {linenum}: {line}"
                    )
                kw, value = line.strip().split(sep=separator3, maxsplit=1)
                result[current_block][kw] = value.strip('"')
            elif any(line.startswith(kw) for kw in allowed_keywords):
                # Ordinary keywords.
                kw, value = line.strip().split(sep=separator3, maxsplit=1)
                result[kw] = value
            else:
                raise ValueError(
                    f"Unexpected line {linenum} in maestro input file: {line}"
                )

        return result
