import os

from pydantic import BaseModel, PrivateAttr
from rdkit import Chem
from copy import deepcopy

from icolos.core.workflow_steps.schrodinger.base import StepSchrodingerBase
from icolos.utils.execute_external.prime import PrimeExecutor
from icolos.utils.execute_external.schrodinger import SchrodingerExecutor

from icolos.core.containers.compound import Conformer

from icolos.utils.enums.program_parameters import PrimeEnum, SchrodingerExecutablesEnum
from icolos.utils.enums.step_enums import StepPrimeEnum, StepGlideEnum
from icolos.core.workflow_steps.step import _LE
from icolos.core.step_utils.sdconvert_util import SDConvertUtil
from icolos.core.step_utils.structcat_util import StructcatUtil
from icolos.utils.general.files_paths import gen_tmp_file
from icolos.utils.general.parallelization import SubtaskContainer, Parallelizer
from tempfile import mkdtemp

_SPE = StepPrimeEnum()
_PE = PrimeEnum()
_SEE = SchrodingerExecutablesEnum()
_SGE = StepGlideEnum()


class StepPrime(StepSchrodingerBase, BaseModel):
    """
    Interface to Schrodinger's Prime mmgbsa implementation
    """

    _schrodinger_executor: SchrodingerExecutor = None

    class Config:
        underscore_attrs_are_private = True

    _sdconvert_util = PrivateAttr()
    _structcat_util = PrivateAttr()

    def __init__(self, **data):
        super().__init__(**data)

        # initialize the executor and test availability
        self._initialize_backend(executor=PrimeExecutor)
        self._check_backend_availability()
        self._schrodinger_executor = SchrodingerExecutor(
            prefix_execution=self.execution.prefix_execution
        )

        # prepare sdconvert utility
        self._sdconvert_util = SDConvertUtil(
            prefix_execution=self.execution.prefix_execution,
            binary_location=self.execution.binary_location,
        )

        # prepare structcat utility
        self._structcat_util = StructcatUtil(
            prefix_execution=self.execution.prefix_execution,
            binary_location=self.execution.binary_location,
        )

    def _execute_prime(self):
        # note, that as the output file name cannot be set (an "-out.maegz" will be attached), this does
        # not need to be heeded here and is encoded in the fixed file name strings

        prime_parallelizer = Parallelizer(func=self._run_subjob)
        n = 1

        while self._subtask_container.done() is False:

            next_batch = self._get_sublists(get_first_n_lists=self._get_number_cores())

            # generate lists for the next batch
            tmp_dirs, complex_paths, output_sdf_paths = self._prepare_batch(next_batch)

            self._apply_token_guard()

            _ = [sub.increment_tries() for element in next_batch for sub in element]
            _ = [sub.set_status_failed() for element in next_batch for sub in element]

            self._logger.log(f"Executing prime for batch {n}", _LE.DEBUG)

            prime_parallelizer.execute_parallel(
                complex_path=complex_paths,
                sdf_output=output_sdf_paths,
                tmp_output_dir=tmp_dirs,
            )

            self._parse_prime_output(
                complex_paths, tmp_dirs, output_sdf_paths, next_batch
            )
            n += 1

    def _parse_prime_output(self, complex_paths, tmp_dirs, output_sdf_paths, batch):
        # go through the batch, get the info from the output file and
        scores = []
        for i in range(len(output_sdf_paths)):
            cur_sublist = batch[i]
            sdf_path = output_sdf_paths[i]
            curr_enum = None
            curr_conformer = None
            mol_supplier = Chem.SDMolSupplier(sdf_path, removeHs=False)
            for mol in mol_supplier:
                # check whether the name corresponds to an enum or conformer
                identifier = str(mol.GetProp("_Name"))
                is_enum = True if len(identifier.split(":")) == 2 else False
                if (
                    not is_enum
                ):  # if we are dealing with a conformer, drop the conformer index to get the enum id
                    enum_index = ":".join(list(mol.GetProp("_Name").split(":"))[:-1])
                else:
                    enum_index = identifier
                # extract the enumeration object, regardless of whether we're dealing with a conformer or an enumeration
                prime_score = mol.GetProp(_SPE.MMGBSA_SCORE)
                for compound in self.get_compounds():
                    for enumeration in compound.get_enumerations():
                        if enumeration.get_index_string() == enum_index:
                            curr_enum = enumeration
                # if we have a conformer, find the conformer with the right ID and append the score to the existing object
                if not is_enum:
                    assert curr_enum is not None
                    for conformer in curr_enum.get_conformers():
                        if conformer.get_index_string() == identifier:
                            curr_conformer = conformer

                else:
                    # if scoring an enumeration, create and attach a conformer out of it
                    curr_conformer = Conformer(conformer=mol)

                    curr_enum.add_conformer(curr_conformer)

                assert curr_conformer is not None
                # now we have the conformer from the originals, set the prime score
                curr_conformer.get_molecule().SetProp(_SPE.MMGBSA_SCORE, prime_score)
                self._logger.log(
                    f"Calculated dG Bind of {prime_score} for conformer {curr_conformer.get_index_string()}",
                    _LE.INFO,
                )
                scores.append(prime_score)

        # after parsing, remove the directories
        self._remove_temporary(tmp_dirs)

        # set success status
        for sublist in batch:
            for task in sublist:
                task.set_status_success()

    def _prepare_batch(self, batch):
        # generate input files for the batch and return tmpdirs

        tmp_dirs = []
        complex_paths = []
        output_sdf_paths = []
        for next_subtask_list in batch:
            tmp_dir = mkdtemp()
            _, tmp_input_sdf_file = gen_tmp_file(suffix=".sdf", dir=tmp_dir)
            _, tmp_input_mae_file = gen_tmp_file(suffix=".maegz", dir=tmp_dir)
            _, tmp_output_sdf_file = gen_tmp_file(suffix=".sdf", dir=tmp_dir)
            writer = Chem.SDWriter(tmp_input_sdf_file)
            for subtask in next_subtask_list:
                mol = deepcopy(subtask.data.get_molecule())
                conf_id = subtask.data.get_index_string()
                mol.SetProp("_Name", conf_id)
                writer.write(mol)
            writer.close()

            # now we have an sdf file with all the conformers from that batch.  Attach the
            structcat_args = [
                "-imae",
                self.settings.additional[_SPE.RECEPTOR],
                "-isd",
                tmp_input_sdf_file,
                "-omae",
                tmp_input_mae_file,
            ]
            self._schrodinger_executor.execute(
                command=_SEE.STRUCTCAT,
                arguments=structcat_args,
                location=tmp_dir,
                check=True,
            )

            tmp_dirs.append(tmp_dir)
            complex_paths.append(tmp_input_mae_file)
            output_sdf_paths.append(tmp_output_sdf_file)

        return tmp_dirs, complex_paths, output_sdf_paths

    def _run_subjob(self, complex_path, sdf_output, tmp_output_dir):

        work_dir = os.getcwd()
        os.chdir(tmp_output_dir)

        arguments = [complex_path, _PE.PRIME_OUTTYPE, _PE.PRIME_OUTTYPE_LIGAND]
        for key in self.settings.arguments.parameters.keys():
            if key not in [_PE.PRIME_OUTTYPE, _PE.PRIME_NJOBS]:
                arguments.append(key)
                arguments.append(str(self.settings.arguments.parameters[key]))
        for flag in self.settings.arguments.flags:
            arguments.append(str(flag))
        if _PE.PRIME_WAIT not in arguments:
            arguments.append(_PE.PRIME_WAIT)

        result = self._backend_executor.execute(
            command=_PE.PRIME_MMGBSA,
            arguments=arguments,
            check=True,
            location=tmp_output_dir,
        )
        print(result)

        output_file = complex_path.split(".")[0] + "-out.maegz"
        assert os.path.isfile(output_file)
        # Convert the mae ligand output back to sdf
        self._sdconvert_util.mae2sdf(output_file, sdf_output)
        os.chdir(work_dir)
        return result

    def execute(self):
        # need to unwrap the conformers to efficiently run in parallel, create lists of subtasks, each with their own files, tmpdirs, then execute them in parallel
        all_conformers = []
        for compound in self.get_compounds():
            for enumeration in compound.get_enumerations():
                if enumeration.get_conformers():
                    # default running mode is to score incoming conformers without changing their configurations
                    for conformer in enumeration.get_conformers():
                        all_conformers.append(conformer)
                else:
                    all_conformers.append(enumeration)

        self._subtask_container = SubtaskContainer(
            max_tries=self.execution.failure_policy.n_tries
        )
        self._subtask_container.load_data(all_conformers)
        self._execute_prime()
        self._logger.log(
            f"Executed Prime for {len(all_conformers)} confomers", _LE.DEBUG
        )
