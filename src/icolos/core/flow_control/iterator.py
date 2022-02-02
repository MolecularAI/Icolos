from typing import Dict, List, Union
from pydantic import BaseModel

from icolos.core.flow_control.flow_control import BaseStepConfig, FlowControlBase
from copy import deepcopy
from icolos.core.job_control.job_control import StepJobControl
from icolos.core.workflow_steps.step import _LE
from icolos.utils.enums.step_enums import StepBaseEnum
from icolos.core.workflow_steps.step import StepBase
from icolos.utils.enums.step_enums import IteratorEnum
import os
from glob import glob

_IE = IteratorEnum
_SBE = StepBaseEnum


class IterSettingsParameters(BaseModel):
    # unpacked version of StepSettingsParameters
    flags: List = []
    parameters: Dict = {}
    additional: Dict = {}
    work_dir: Union[List, str] = []


class IterParallelizer(BaseModel):

    # if true, steps must be totally independent, the iterator step
    parallelize: bool = False
    cores: int = 1
    dependent_steps: int = None


class IterSettings(BaseModel):
    # dictionary of settings to change

    # settings: {step_id: {IterSettingsParameters}}
    settings: Dict[str, IterSettingsParameters] = {}
    iter_mode: _IE = _IE.N_ITERS
    n_iters: int = None
    parallelizer_settings: IterParallelizer = IterParallelizer()


class StepIterator(FlowControlBase, BaseModel):
    """
    Implements iterator mechanism:
    wraps one or multiple steps and generates n copies of that set of steps according to iter_config
    Becomes master job when parallize=True, using  Icolos JobControl to interface with external resources
    """

    # holds the dict of iterables for the bits to chang
    iter_settings: IterSettings = IterSettings()

    def __init__(self, **data):
        super().__init__(**data)
        # when init_step_from_dict calls this method, we need to initialise a list of steps,
        # controlled by iter_settings.iter_mode
        self.initialized_steps = self._initialize_steps()
        # either generate a list, if serial execution, or initialize a single JobControl
        # step with each config as an initialized step

    def _initialize_n_iters(self) -> List:
        """
        Initialise n identical copies of the same step config
        """
        init_steps = []
        for i in range(self.iter_settings.n_iters):

            list_step_conf = deepcopy(self.base_config)

            # hand all steps over to the config updater
            formatted_confs = self._update_config(list_step_conf, f"run_{i}")
            for step_conf in formatted_confs:
                initialized_step = self._initialize_step_from_dict(step_conf._as_dict())
                init_steps.append(initialized_step)
        return init_steps

    def _modify_settings(
        self, settings: BaseStepConfig, step_config: BaseStepConfig, i: int
    ):
        base_conf = deepcopy(step_config)
        iter_settings = deepcopy(settings)
        if settings.flags:
            settings.flags = {"flags": settings.flags}

            # iterate over the flags
            # if it's been converted, hence there are flags to be applied
            base_conf.settings.arguments.flags.append(iter_settings.flags.values()[i])
        for (
            key,
            val,
        ) in iter_settings.parameters.items():
            base_conf.settings.arguments.parameters[key] = val[i]
        for (
            key,
            val,
        ) in iter_settings.additional.items():
            # however many lists of n items
            base_conf.settings.additional[key] = val[i]

        # in pmx abfe, we iterate over workdirs
        if isinstance(iter_settings.work_dir, str):
            # a single top dir has been provided, extract the subdirs and walk over these
            work_dirs = glob(f"{iter_settings.work_dir}/**/", recursive=True)
            # TODO, automatically acquire the number of iters from this
            base_conf.work_dir = os.path.join(iter_settings.work_dir, work_dirs[i])
        elif isinstance(iter_settings.work_dir, list) and iter_settings.work_dir:
            base_conf.work_dir = iter_settings.work_dir[i]

        return base_conf

    def _initialize_settings(self) -> List:
        """
        Iterate through all settings step-wise, changing all setting blocks simultaneously, returning n initialised steps for n
        """
        init_steps = []
        for i in range(self.iter_settings.n_iters):

            # iterate over the steps in the base config, and the corresponding settings, if these are to be modified
            step_sublist = []
            for step_config in self.base_config:

                # check if we need to iterate through settings in this step, else just use the base config
                if step_config.step_id in self.iter_settings.settings.keys():
                    settings = self.iter_settings.settings[step_config.step_id]
                    step_sublist.append(self._modify_settings(settings, step_config, i))
                else:
                    step_sublist.append(step_config)

            # update all configs with references to updated step_ids etc
            formatted_configs = self._update_config(step_sublist, f"run_{i}")
            for step_conf in formatted_configs:
                initialized_step = self._initialize_step_from_dict(step_conf._as_dict())
                init_steps.append(initialized_step)
        return init_steps

    # def _initialize_compounds(self):
    #     """
    #     Generates n copies of a step, each with a single compound loaded from the source step
    #     * Only the first step in base_config needs updating, downstream data handover from this step is handed properly anyway
    #     """
    #     init_steps = []
    #     # TODO: get the number of compounds automatically?
    #     for i in range(self.iter_settings.n_iters):
    #         list_step_conf = deepcopy(self.base_config)
    #         first_step_config = list_step_conf[0]
    #         # probably only expecting one set of input compounds but this will select the ith for all inputs
    #         for inp_block in first_step_config.input.compounds:
    #             inp_block.selected_compound_id = i
    #         formatted_confs = self._update_config(list_step_conf, f"run_{i}")
    #         for step_conf in formatted_confs:
    #             initialized_step = self._initialize_step_from_dict(step_conf.as_dict())
    #             init_steps.append(initialized_step)
    #     return init_steps

    def _initialize_steps(self) -> Union[List, StepBase]:
        """
        Handle step init according to config
        Returns a list of steps if serial execution (default)
        Returns a Step-like JobControl object if parallelisation is specified.
        """
        steps = []
        if self.iter_settings.iter_mode == _IE.N_ITERS:
            # simplest mode, just n repeats of the same step
            steps += self._initialize_n_iters()

        elif self.iter_settings.iter_mode == _IE.SINGLE:
            # for n different settings, iterate through each, returning n steps
            steps += self._initialize_settings()
        elif self.iter_settings.iter_mode == _IE.ALL:
            raise NotImplementedError
            # initialise all combinations of steps by combining settings
            # steps.append(self._initialize_combined())

        self._logger.log(
            f"Iterator has initialized {len(steps)} steps for step  {self.base_config[0].step_id}",
            _LE.DEBUG,
        )
        if not self.iter_settings.parallelizer_settings.parallelize:
            return steps
        else:

            wrapper = StepJobControl(
                step_id="JobControl",
                type=_SBE.STEP_JOB_CONTROL,
                initialized_steps=steps,
                parallel_execution=self.iter_settings.parallelizer_settings,
            )
            return wrapper

    def _update_config(
        self, step_conf: List[BaseStepConfig], run_id: str
    ) -> List[BaseStepConfig]:
        """
        Manages modifications to each step in the config:
        * step_id is updated with the run_id
        * any references to other step_ids (e.g. in input) contained in the base config are updated to reflect the change
        * writeout paths are updated to separate output from each of the runs
        """
        original_step_ids = [conf.step_id for conf in step_conf]
        formatted_confs = []
        for conf in step_conf:
            # modify the step_id
            st_id = conf.step_id
            conf.step_id = st_id + "_" + run_id
            # modify the writeout paths: add a key_value dir the writeout path
            for idx, block in enumerate(conf.writeout):
                if block.destination.resource is not None:
                    resource = block.destination.resource
                    parts = resource.split("/")
                    new_resource = os.path.join("/".join(parts[:-1]), run_id, parts[-1])
                    block.destination.resource = new_resource

            # modify the step_input blocks if they reference a step_id contained in step_conf
            # treat compounds
            for comp in conf.input.compounds:
                if comp.source in original_step_ids:
                    comp.source += f"_{run_id}"

            for gen in conf.input.generic:
                if gen.source in original_step_ids:
                    gen.source += f"_{run_id}"

            # TODO: this is a bodge for now
            # we have an edge case in data manipularion that needs to match compounds those from another step, the source name needs the same treatment
            if conf.type.upper() == _SBE.STEP_DATA_MANIPULATION:
                if (
                    _SBE.INPUT_SOURCE in conf.settings.additional.keys()
                    and conf.settings.additional[_SBE.INPUT_SOURCE] in original_step_ids
                ):

                    conf.settings.additional[_SBE.INPUT_SOURCE] += f"_{run_id}"

            formatted_confs.append(conf)

        return formatted_confs
