class ExecutionFailed(Exception):
    pass


class StepFailed(Exception):
    pass


class ContainerCorrupted(Exception):
    pass


def get_exception_message(e: Exception):
    if e is None:
        return None
    if hasattr(e, "message"):
        return e.message
    else:
        return e


def get_exception_type(e: Exception) -> str:
    if e is None:
        return None
    return type(e).__name__
