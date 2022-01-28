from copy import deepcopy
from typing import List, Tuple
from icolos.core.workflow_steps.step import StepBase
from pydantic import BaseModel
import tempfile
from icolos.utils.enums.step_enums import StepCressetEnum
from icolos.utils.execute_external.cresset_executor import CressetExecutor
from icolos.utils.general.files_paths import gen_tmp_file
from icolos.core.workflow_steps.step import _LE
import os
from icolos.utils.general.parallelization import Parallelizer, SubtaskContainer


_SCE = StepCressetEnum()


class StepCressetEC(StepBase, BaseModel):
    def __init__(self, **data):
        super().__init__(**data)

        self._initialize_backend(executor=CressetExecutor)
        self._check_backend_availability()

    def _prepare_tmp_input(self, batch: List) -> Tuple[List, List]:
        conformers = []
        tmp_dirs = []
        protein = self.data.generic.get_argument_by_extension(
            "pdb", rtn_file_object=True
        )
        for sublist in batch:
            for task in sublist:
                conformer = task.data
                conformers.append(conformer)

                # generate the tmpdir
                tmp_dir = tempfile.mkdtemp()
                tmp_dirs.append(tmp_dir)
                _, path_input_sdf = gen_tmp_file(
                    prefix="tmp_", suffix=".sdf", dir=tmp_dir
                )
                conformer.write(path=path_input_sdf)

                # write the protein to that tmpdir
                protein.write(path=os.path.join(tmp_dir, "protein.pdb"), join=False)

        return conformers, tmp_dirs

    def _execute_cresset_ec_parallel(self):
        parallelizer = Parallelizer(func=self._run_conformer)
        n = 1

        while self._subtask_container.done() is False:

            next_batch = self._get_sublists(
                get_first_n_lists=self._get_number_cores()
            )  # return n lists of length max_sublist_length
            _ = [sub.increment_tries() for element in next_batch for sub in element]
            _ = [sub.set_status_failed() for element in next_batch for sub in element]

            conformers, tmp_dirs = self._prepare_tmp_input(next_batch)
            self._logger.log(
                f"Executing Cresset EC for batch {n} containing {len(conformers)} conformers",
                _LE.DEBUG,
            )

            parallelizer.execute_parallel(tmp_dir=tmp_dirs, conformer=conformers)

            results = self._parse_results(tmp_dirs, conformers)

            for sublist, result in zip(next_batch, results):
                # TODO: this only works if max length sublist == 1, fine for now as that is all turbomole can handle
                for task in sublist:
                    if result == _SCE.SUCCESS:
                        task.set_status_success()
                    else:
                        task.set_status_failed()
            self._remove_temporary(tmp_dirs)
            n += 1

    def _parse_results(self, tmp_dirs: List, conformers: List):
        # walk over the directory structure, parse the output file, identify the conformer, attach a tag to the mol object
        # TODO: No idea what the output looks like for this, write the parser!!
        pass

    def execute(self):
        # unroll all conformers
        all_conformers = []
        for compound in self.get_compounds():
            for enum in compound.get_enumerations():
                if self._input_object_empty(enum):
                    continue
                else:
                    for conformer in enum.get_conformers():
                        conf = deepcopy(conformer)
                        all_conformers.append(conf)

        self._subtask_container = SubtaskContainer(
            max_tries=self.execution.failure_policy.n_tries
        )
        self._subtask_container.load_data(all_conformers)
        self._execute_cresset_ec_parallel()

    def _run_conformer(self):
        # run a single conformer through Flare's EC
        self._backend_executor.execute()

    # execution is
    # module load Flare && pyflare electrostaticcomplementarity.py -p protein.pdb ligands.sdf
