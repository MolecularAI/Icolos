from icolos.utils.enums.program_parameters import GromacsEnum
from icolos.utils.execute_external.execute import ExecutorBase

_GE = GromacsEnum()


class GromacsExecutor(ExecutorBase):
    def __init__(self, prefix_execution=None, binary_location=None):
        super().__init__(
            prefix_execution=prefix_execution, binary_location=binary_location
        )

    def execute(
        self, command: str, arguments: list, check=True, location=None, pipe_input=None
    ):
        gmx_commands = [
            _GE.EDITCONF,
            _GE.GENION,
            _GE.GROMPP,
            _GE.SOLVATE,
            _GE.MDRUN,
            _GE.MPI_MDRUN,
            _GE.PDB2GMX,
            _GE.MAKE_NDX,
            _GE.GENRESTR,
            _GE.TRJCONV,
            _GE.TRJCAT,
            _GE.CLUSTER,
            _GE.MMPBSA,
            _GE.DO_DSSP,
            _GE.RMS,
        ]

        if not any([cmd in command for cmd in gmx_commands]):
            raise ValueError(
                "Command must be present in internal list of GROMACS executables"
            )
        # handle for dealing with programs that want interactive input
        return super().execute(
            command=command,
            arguments=arguments,
            check=check,
            location=location,
            pipe_input=pipe_input,
        )

    def is_available(self):
        try:
            result = self.execute(command=_GE.PDB2GMX, arguments=[], check=False)
            if _GE.PDB2GMX_FAIL_ID_STRING in result.stderr:
                return True
            return False
        except Exception as e:
            return False
