from pydantic import BaseModel
from icolos.core.containers.gmx_state import GromacsState
from icolos.core.workflow_steps.step import StepBase
from typing import List
from icolos.utils.enums.step_enums import StepCavExploreEnum

_SFP = StepCavExploreEnum()


class StepFpocketBase(StepBase, BaseModel):
    eps: float = None
    iso_value: int = None
    threshold: float = None
    min_samples: int = None
    format_: str = None

    def __init__(self, **data):
        super().__init__(**data)

    def _write_input_files(self, tmp_dir, topol: GromacsState):
        for file in self.data.generic.get_flattened_files():
            file.write(tmp_dir)
        if topol is not None:
            # handle structure files
            if not self.data.generic.get_files_by_extension("gro"):
                if topol.structures:
                    topol.write_structure(tmp_dir)
            if not self.data.generic.get_files_by_extension("top"):
                if topol.top_lines:
                    topol.write_topol(tmp_dir)
            if not self.data.generic.get_files_by_extension("tpr"):
                if topol.tprs:
                    topol.write_tpr(tmp_dir)
            if not self.data.generic.get_files_by_extension("xtc"):
                if topol.trajectories:
                    topol.write_trajectory(tmp_dir)
            if not self.data.generic.get_files_by_extension("ndx"):
                if topol.ndx:
                    topol.write_ndx(tmp_dir)

    def _parse_arguments(self, flag_dict: dict, args: list = None) -> List:
        arguments = args if args is not None else []
        # first add the settings from the command line
        for key in self.settings.arguments.parameters.keys():
            arguments.append(key)
            arguments.append(str(self.settings.arguments.parameters[key]))
        for flag in self.settings.arguments.flags:
            arguments.append(str(flag))
        for key, value in flag_dict.items():
            # only add defaults if they have not been specified in the json
            if key not in arguments:
                arguments.append(key)
                arguments.append(value)
        return arguments

    def _set_mdpocket_args(self):
        if self.settings.additional is not None:
            keys = self.settings.additional.keys()

            self.eps = self.settings.additional[_SFP.EPS] if _SFP.EPS in keys else 3
            self.iso_value = (
                self.settings.additional[_SFP.ISO_VALUE]
                if _SFP.ISO_VALUE in keys
                else 0.5
            )
            self.threshold = (
                self.settings.additional[_SFP.THRESHOLD]
                if _SFP.THRESHOLD in keys
                else 20.0
            )
            self.min_samples = (
                self.settings.additional[_SFP.MIN_SAMPLES]
                if _SFP.MIN_SAMPLES in keys
                else 25
            )
            if _SFP.TRAJ_TYPE in keys:
                if self.settings.additional[_SFP.TRAJ_TYPE].lower() == "gromacs":
                    self.format_ = "xtc"
                elif self.settings.additional[_SFP.TRAJ_TYPE].lower() == "desmond":
                    self.format_ = "dtr"
                else:
                    raise ValueError(
                        "Only Desmond and GROMACS trajectory types are supported"
                    )
            else:
                raise ValueError("Trajectory format was not specified!")
