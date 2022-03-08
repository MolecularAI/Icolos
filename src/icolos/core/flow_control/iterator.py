from typing import Dict, List, Union
from pydantic import BaseModel
from icolos.core.composite_agents.workflow import WorkFlow

from icolos.core.flow_control.flow_control import BaseStepConfig, FlowControlBase
from copy import deepcopy
from icolos.core.workflow_steps.step import _LE
from icolos.core.step_dispatch.dispatcher import StepDispatcher
from icolos.utils.enums.composite_agents_enums import WorkflowEnum
from icolos.utils.enums.step_enums import StepBaseEnum
from icolos.core.workflow_steps.step import StepBase
from icolos.utils.enums.step_enums import IteratorEnum
from icolos.core.composite_agents.workflow import WorkflowHeaderParameters
import os

_IE = IteratorEnum
_SBE = StepBaseEnum
_WE = WorkflowEnum()


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
    remove_temporary_files: bool = True
    single_directory: bool = False
    parallelizer_settings: IterParallelizer = IterParallelizer()


class StepIterator(FlowControlBase, BaseModel):
    """
    Implements iterator mechanism:
    wraps one or multiple steps and generates n copies of that set of steps according to iter_config
    Becomes master job when parallize=True, using  Icolos JobControl to interface with external resources
    """

    # holds the dict of iterables for the bits to chang
    iter_settings: IterSettings = IterSettings()
    dispatcher: StepDispatcher = None

    def __init__(self, **data):
        super().__init__(**data)

        self.dispatcher = self._initialize_workflows()

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

        return base_conf

    def _initialize_n_iters(self) -> List[WorkFlow]:
        """
        Initialise n identical copies of the same step config
        """
        workflows = []
        for i in range(self.iter_settings.n_iters):
            workflow_steps = []
            list_step_conf = deepcopy(self.base_config)

            # hand all steps over to the config updater
            step_confs = self._update_config(list_step_conf, f"run_{i}")
            for step in step_confs:
                workflow_steps.append(step.as_dict())
            wf_config = {
                _WE.HEADER: {
                    _WE.ID: f"workflow_{i}",
                    _WE.ENVIRONMENT: {_WE.ENVIRONMENT_EXPORT: []},
                    _WE.GLOBAL_VARIABLES: {},
                    _WE.GLOBAL_SETTINGS: {
                        "remove_temporary_files": self.iter_settings.remove_temporary_files,
                        "single_directory": self.iter_settings.single_directory,
                    },
                },
                _WE.STEPS: workflow_steps,
            }
            # manually set some things for these workflows

            workflows.append(WorkFlow(**wf_config))
        return workflows

    def _initialize_settings(self) -> List:
        """
        Iterate through all settings step-wise, changing all setting blocks simultaneously, returning n initialised steps for n
        """
        init_steps = []
        for i in range(self.iter_settings.n_iters):

            # iterate over the steps in the base config, and the corresponding settings, if these are to be modified
            step_sublist = []
            for step_config in deepcopy(self.base_config):

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

    def _initialize_workflows(self) -> StepDispatcher:
        """
        Handle step init according to config
        Returns a Step-like JobControl that holds a list of separate workflows, which will be executed according to the parallelization scheme.
        """
        workflows = []
        if self.iter_settings.iter_mode == _IE.N_ITERS:
            # simplest mode, generate n independent workflows with the same set of steps
            workflows += self._initialize_n_iters()

        elif self.iter_settings.iter_mode == _IE.SINGLE:
            raise NotImplementedError
            # for n different settings, iterate through each, returning n steps
            steps += self._initialize_settings()
        else:
            raise NotImplementedError

        self._logger.log(
            f"Initialized {len(workflows)} jobs for step  {self.base_config[0].step_id}",
            _LE.DEBUG,
        )
        dispatcher = StepDispatcher(
            step_id="step_dispatcher",
            type=_SBE.STEP_DISPATCHER,
            workflows=workflows,
            parallel_execution=self.iter_settings.parallelizer_settings,
        )
        return dispatcher

    def _update_config(
        self, step_conf: List[BaseStepConfig], run_id: str
    ) -> List[BaseStepConfig]:
        """
        Manages modifications to each step in the config:
        * writeout paths are updated to separate output from each of the runs
        """

        formatted_confs = []
        for conf in step_conf:

            # modify the writeout paths: add a key_value dir the writeout path
            for block in conf.writeout:
                if block.destination.resource is not None:
                    resource = block.destination.resource
                    parts = resource.split("/")
                    new_resource = os.path.join("/".join(parts[:-1]), run_id, parts[-1])
                    block.destination.resource = new_resource

            formatted_confs.append(conf)

        return formatted_confs
