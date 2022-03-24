from typing import List
from pydantic.main import BaseModel
from icolos.core.composite_agents.workflow import WorkFlow
from icolos.core.workflow_steps.step import StepBase
from icolos.utils.general.parallelization import Parallelizer, SubtaskContainer
from icolos.core.workflow_steps.step import _LE


class IterParallelizer(BaseModel):
    # config block controlling how the steps are parallelized
    # if you are executing a 5 step workflow with 10 repeats, dependent_steps = 5, cores = 10
    # this will allow each independent replica to be allocated to a single job queue, retaining step order
    parallelize: bool = False
    jobs: int = 1
    dependent_steps: int = None


class StepDispatcher(StepBase, BaseModel):
    """
    Step class containing job control functionality required for StepIterator, supports Slurm for job scheduling
    Supports running Icolos process as master job for parallel step execution on cluster.  Generates a pool of initialized steps to be executed, based on the

    Step-type class for disaptching multiple steps in parallel, useful for executing multiple batch jobs simultaneously
    """

    workflows: List = []
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
        # Spin up multiple processes.
        self.execution.parallelization.jobs = self.parallel_execution.jobs

        # TODO, we can repeat entire workflows if we want, I'm not sure this makes sense though
        self._subtask_container = SubtaskContainer(max_tries=1)
        self._subtask_container.load_data(self.workflows)

        parallelizer = Parallelizer(func=self.execute_workflow)
        n = 1

        while self._subtask_container.done() is False:

            next_batch = self._get_sublists(
                get_first_n_lists=self.parallel_execution.jobs
            )  # return n lists of length max_sublist_length
            _ = [sub.increment_tries() for element in next_batch for sub in element]
            _ = [sub.set_status_failed() for element in next_batch for sub in element]

            self._logger.log(
                f"Starting {len(next_batch)} parallel jobs under Icolos JobControl, execution batch {n}",
                _LE.INFO,
            )

            jobs = self._prepare_batch(next_batch)

            result = parallelizer.execute_parallel(jobs=jobs)

            # TODO: sucessful execution of each step is not explicitly checked,
            # the step is responsible for throwing errors if something has gone wrong
            for task in next_batch:
                for subtask in task:
                    subtask.set_status_success()

    def execute_workflow(self, jobs):
        # submits then monitors the step
        wf_data = self.get_workflow_object().workflow_data
        for idx, job in enumerate(jobs):
            # copy existing wf data up to this point int othe new wf object
            job.initialize()
            job.workflow_data = wf_data
            self._logger.log(f"Executing workflow {idx} of {len(jobs)}", _LE.DEBUG)
            job.execute()
