from icolos.core.workflow_steps.step import StepBase
from pydantic import BaseModel
from icolos.utils.execute_external.rosetta import RosettaExecutor

# Wrapping for Rosetta ab initio structure prediction
# Note the execution is finicky, and requires some proper set up to get this to work
# we run everything locally, since public web servers are not an option for us.
# some config is required to get make_fragments.pl to run with all its dependencies in place


class StepRosettaAbinitio(StepBase, BaseModel):
    def __init__(self, **data):
        super().__init__(data)

        self._inititalize_backend(executor=RosettaExecutor)
