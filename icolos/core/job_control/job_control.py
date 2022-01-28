from typing import List
from pydantic.main import BaseModel

from icolos.core.workflow_steps.step import StepBase
from icolos.utils.general.parallelization import Parallelizer, SubtaskContainer
from icolos.core.workflow_steps.step import _LE


class IterParallelizer(BaseModel):
    # config block controlling how the steps are parallelized
    # if you are executing a 5 step workflow with 10 repeats, dependent_steps = 5, cores = 10
    # this will allow each independent replica to be allocated to a single job queue, retaining step order
    parallelize: bool = False
    cores: int = 1
    dependent_steps: int = None


class StepJobControl(StepBase, BaseModel):
    """
    Step class containing job control functionality required for StepIterator, supports Slurm for job scheduling
    Supports running Icolos process as master job for parallel step execution on cluster.  Generates a pool of initialized steps to be executed, based on the
    """

    initialized_steps: List = []
    # expect the parallel execution block to be handed over from flow control
    parallel_execution: IterParallelizer = IterParallelizer()

    def __init__(self, **data):
        super().__init__(**data)

    def _prepare_batch(self, batch) -> List[List[StepBase]]:

        batch_steps = []
        for sublist in batch:
            sublist_steps = []
            for task in sublist:
                sublist_steps.append(task.data)
            batch_steps.append(sublist_steps)
        return batch_steps

    def execute(self):
        """
        Execute multiple steps in parallel
        """
        # Spin up multiple processes
        self.execution.parallelization.cores = self.parallel_execution.cores
        # each subtask needs to contain an entire mini workflow to be executed sequentially,
        self.execution.parallelization.max_length_sublists = (
            self.parallel_execution.dependent_steps
        )

        # if we try steps multiple times, we have steps fail depending on its dependency on a
        # previous step - too complicated
        self._subtask_container = SubtaskContainer(max_tries=1)
        self._subtask_container.load_data(self.initialized_steps)

        parallelizer = Parallelizer(func=self._run_step)
        n = 1

        while self._subtask_container.done() is False:

            next_batch = self._get_sublists(
                get_first_n_lists=self.parallel_execution.cores
            )  # return n lists of length max_sublist_length
            _ = [sub.increment_tries() for element in next_batch for sub in element]
            _ = [sub.set_status_failed() for element in next_batch for sub in element]

            self._logger.log(
                f"Starting {len(next_batch)} parallel jobs under Icolos JobControl, execution batch {n}",
                _LE.INFO,
            )

            steps = self._prepare_batch(next_batch)

            result = parallelizer.execute_parallel(steps=steps)

            # sucessful execution of each step is not explicitly checked,
            # the step is responsible for throwing errors if something has gone wrong
            for task in next_batch:
                for subtask in task:
                    subtask.set_status_success()

    def _run_step(self, steps: List[StepBase]):
        # submits then monitors the step
        for step in steps:  # length max_len_sublist
            # at this point the internal steps don't have their data initialised
            step.generate_input()
            step.execute()
            step.process_write_out()
