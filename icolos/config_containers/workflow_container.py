from icolos.config_containers.container import ConfContainer


class WorkflowContainer(ConfContainer):
    def __init__(self, conf, validation=True):
        super().__init__(conf=conf)

        # TODO: include validation with JSON Schema
        if validation:
            self.validate()

    def validate(self):
        pass
