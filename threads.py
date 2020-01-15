import logging
import time
from abc import (
    ABCMeta,
    abstractmethod,
)
from concurrent.futures import (
    Future,
    ThreadPoolExecutor,
    as_completed,
)
from typing import (
    Iterable,
    List,
    Optional,
)

logger = logging.getLogger(__name__)


class DeferredTaskExecutor:

    def __init__(self, max_workers: Optional[int] = None) -> None:
        super().__init__()
        self.tpe = ThreadPoolExecutor(max_workers=max_workers)
        self.futures = set()
        self._errors = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._errors = self._collect_futures()
        self.tpe.shutdown(wait=True)
        return False

    @property
    def errors(self):
        assert self._errors is not None, 'Errors can only be used after context manager has been exited.'
        return self._errors

    def defer(self,
              callable_,
              *args,
              run_after: Optional[Iterable[Future]] = None,
              start_time: Optional[float] = None,
              delay: Optional[float] = None,
              **kwargs) -> Future:
        """
        Invoke the given callable (typically a method of this class or a function nested in a method) with the given
        arguments and after the preconditions are met.

        :param callable_: the callable to invoke

        :param args: the positional arguments to pass to the callable

        :param kwargs: the keyword arguments to pass to the callable

        :param run_after: the futures representing other callables that must complete successfully before
                          this callable is invoked

        :param start_time: an optional absolute point in time (as returned by time.time())
                           before which that task will not be invoked, defaults to now

        :param delay: an optional number of seconds that will be added to start_time

        :return: a Future instance representing the callable
        """
        if start_time is None:
            if delay is not None:
                start_time = time.time() + delay
        else:
            if delay is not None:
                start_time = start_time + delay

        def run_if_possible():
            can_run = self._check_run_after(run_after) if run_after else True
            if can_run is False:
                raise self.UnsatisfiedDependency()
            elif can_run is True and (start_time is None or start_time < time.time()):
                return callable_(*args, **kwargs)
            else:
                return self.defer(callable_, *args, run_after=run_after, start_time=start_time, **kwargs)

        def log_exceptions_early(future):
            e = future.exception()
            if e is not None and not isinstance(e, self.UnsatisfiedDependency):
                logger.warning('Exception in deferred callable', exc_info=True)

        future = self.tpe.submit(run_if_possible)
        future.add_done_callback(log_exceptions_early)
        self.futures.add(future)
        return future

    class UnsatisfiedDependency(RuntimeError):
        pass

    def _check_run_after(self, run_after: Iterable[Future]) -> Optional[bool]:
        for future in run_after:
            while True:
                if future.done():
                    if future.exception():
                        return False  # at least one future failed
                    else:
                        # Tasks that call _defer() will return a future which needs to be examined recursively.
                        # This tail recursion could be of arbitrary depth.
                        result = future.result()
                        if isinstance(result, Future):
                            future = result
                        else:
                            break
                else:
                    return None  # some futures are not yet done
        return True  # all futures succeeded

    def _collect_futures(self) -> List[BaseException]:
        errors = []
        num_secondary_errors = 0
        while self.futures:
            for future in as_completed(self.futures):
                e = future.exception()
                if e is not None:
                    if isinstance(e, self.UnsatisfiedDependency):
                        num_secondary_errors += 1
                    else:
                        errors.append(e)
                self.futures.remove(future)
        # We cannot have any secondary errors without primary ones
        assert bool(errors) or not bool(num_secondary_errors)
        return errors
