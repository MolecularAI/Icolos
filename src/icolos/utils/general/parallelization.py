import math
import multiprocessing
from typing import List, Callable, Dict, Any
from pydantic import BaseModel
from icolos.utils.enums.parallelization import ParallelizationEnum

_PE = ParallelizationEnum


class Subtask(BaseModel):
    status: _PE = _PE.STATUS_READY
    times_tried: int = 0
    data: Any

    def increment_tries(self):
        self.times_tried += 1

    def set_status(self, status: str):
        self.status = status

    def set_status_failed(self):
        self.set_status(_PE.STATUS_FAILED)

    def set_status_success(self):
        self.set_status(_PE.STATUS_SUCCESS)


class SubtaskContainer(BaseModel):
    max_tries: int
    subtasks: List[Subtask] = []

    def __init__(self, **data):
        super().__init__(**data)

    def clear(self):
        self.subtasks = []

    def load_data(self, data: List[Any]):
        self.clear()
        self.add_data(data=data)

    def add_data(self, data: List[Any]):
        for data_element in data:
            self.subtasks.append(
                Subtask(status=_PE.STATUS_READY, times_tried=0, data=data_element)
            )

    def get_todo_tasks(self) -> List[Subtask]:
        todo_subtasks = []
        for subtask in self.subtasks:
            if (
                subtask.status == _PE.STATUS_READY
                or subtask.status == _PE.STATUS_FAILED
            ) and subtask.times_tried < self.max_tries:
                todo_subtasks.append(subtask)
        return todo_subtasks

    def get_done_tasks(self) -> List[Subtask]:
        done_subtasks = []
        for subtask in self.subtasks:
            if (
                subtask.status == _PE.STATUS_SUCCESS
                or subtask.times_tried >= self.max_tries
            ):
                done_subtasks.append(subtask)
        return done_subtasks

    def get_sublists(
        self, partitions=None, slice_size=None, get_first_n_lists=None
    ) -> List[List[Subtask]]:
        if partitions is None and slice_size is None:
            raise ValueError("Either specify partitions or slice size.")

        # only get tasks that are not yet completed or have some tries left
        subtasks = self.get_todo_tasks()

        # decide on the chunk size, either by doing partitions or by specifying the slice size directly
        sublists = []
        if partitions is not None:
            chunk_size = int(math.ceil(len(subtasks) / partitions))
        else:
            chunk_size = slice_size

        # wrap the tasks in lists as required
        for i in range(0, len(subtasks), chunk_size):
            sublist = []
            for j in range(i, min(i + chunk_size, len(subtasks))):
                sublist.append(subtasks[j])
            sublists.append(sublist)

        if get_first_n_lists is not None and len(sublists) > get_first_n_lists:
            return sublists[:get_first_n_lists]
        else:
            return sublists

    def done(self) -> bool:
        for subtask in self.subtasks:
            if subtask.status == _PE.STATUS_SUCCESS:
                continue
            if subtask.status == _PE.STATUS_READY or (
                subtask.status == _PE.STATUS_FAILED
                and subtask.times_tried < self.max_tries
            ):
                return False
        return True

    def any_failed(self) -> bool:
        if len(
            [True for subtask in self.subtasks if subtask.status == _PE.STATUS_FAILED]
        ):
            return True
        return False

    def set_max_tries(self, max_tries: int):
        self.max_tries = max_tries

    def __len__(self) -> int:
        return len(self.subtasks)


class Parallelizer(BaseModel):
    func: Callable

    def __init__(self, **data):
        super().__init__(**data)

    def rearrange_input(self, inp_dict: Dict[str, List]) -> List[Dict]:
        return [dict(zip(inp_dict, ele)) for ele in zip(*inp_dict.values())]

    def execute_parallel(self, **kwargs):
        # translate the dictionary with the lists of arguments into a list of individual dictionaries
        # e.g. {'one': [1, 2, 3], 'two': ['aaaa', 'bbb', 'cc'], 'three': [0.2, 0.2, 0.1]} --->
        # [{'one': 1, 'two': 'aaaa', 'three': 0.2},
        #  {'one': 2, 'two': 'bbb', 'three': 0.2},
        #  {'one': 3, 'two': 'cc', 'three': 0.1}]
        list_exec = self.rearrange_input(kwargs)

        # # run in parallel; wait for all subjobs to finish before proceeding
        processes = []
        for subprocess_args in list_exec:
            p = multiprocessing.Process(target=self.func, kwargs=subprocess_args)
            processes.append(p)
            p.start()

        for p in processes:
            p.join()
