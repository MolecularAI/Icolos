import functools
import time
from typing import Any
from pydantic import BaseModel


class RetryResult(BaseModel):
    success: bool
    tries: int
    result: Any = None
    exception: Exception = None

    class Config:
        arbitrary_types_allowed = True

    def __init__(self, **data):
        super().__init__(**data)


# TODO: do a unit test for this
def retry(n_tries, retry_wait_seconds, allowed_exceptions=()):
    if n_tries < 1:
        n_tries = 1

    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs) -> RetryResult:
            exc = None
            for idx in range(n_tries):
                try:
                    result = func(*args, **kwargs)
                    return RetryResult(
                        success=True, tries=idx + 1, result=result, exception=None
                    )
                except allowed_exceptions as e:
                    exc = e
                time.sleep(retry_wait_seconds)
            return RetryResult(success=False, tries=n_tries, result=None, exception=exc)

        return wrapper

    return decorator
