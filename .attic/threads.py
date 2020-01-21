from collections import defaultdict
from concurrent.futures import (
    Future,
    ThreadPoolExecutor,
    as_completed,
)
import logging
from random import random
from threading import Lock
import time
from typing import (
    Any,
    Callable,
    Hashable,
    Iterable,
    List,
    MutableMapping,
    Optional,
)

logger = logging.getLogger(__name__)

Key = Hashable
ErrorsByKey = MutableMapping[Key, List[BaseException]]


class DeferredFuture:

    def __init__(self,
                 executor: 'DeferredTaskExecutor',
                 key: Key,
                 parents: Optional[Iterable['DeferredFuture']],
                 *args: Any,
                 **kwargs: Any) -> None:
        super().__init__()
        self.children: Optional[List['DeferredFuture']] = []
        self.child_lock = Lock()
        self.num_parents = 0
        self.executor = executor
        self.key = key
        self.args = args
        self.kwargs = kwargs
        if parents is None:
            self._schedule()
        else:
            parents = list(parents)
            if parents:
                self.num_parents = len(parents)
                for parent in parents:
                    parent._add_child(self)
            else:
                self._schedule()

    def _add_child(self, child: 'DeferredFuture'):
        with self.child_lock:
            if self.children is None:
                # already done, notify immediately
                child._notify()
            else:
                self.children.append(child)

    def _schedule(self):
        # noinspection PyProtectedMember
        future = self.executor._submit(self)
        future.add_done_callback(self._done_callback)

    def _notify(self):
        assert self.num_parents > 0
        self.num_parents -= 1
        if self.num_parents == 0:
            self._schedule()

    def _done_callback(self, future):
        e = future.exception()
        if e is None:
            with self.child_lock:
                for child in self.children:
                    # noinspection PyProtectedMember
                    child._notify()
                # signal that we're done
                self.children = None
        else:
            logger.warning('Exception in deferred callable', exc_info=e)


class DeferredTaskExecutor:

    def __init__(self, max_workers: Optional[int] = None) -> None:
        super().__init__()
        self._tpe = ThreadPoolExecutor(max_workers=max_workers)
        self._futures: MutableMapping[Future, Key] = {}
        self._errors: Optional[ErrorsByKey] = None

    def __enter__(self):
        return self

    def defer(self,
              key: Key,
              callable_: Callable,
              *args,
              run_after: Optional[Iterable[DeferredFuture]] = None,
              **kwargs) -> DeferredFuture:
        """
        Invoke the given callable (typically a method of this class or a function nested in a method) with the given
        arguments and after the preconditions are met.

        :param key: a key by which to identify this callable when reporting exceptions

        :param callable_: the callable to invoke

        :param args: the positional arguments to pass to the callable

        :param kwargs: the keyword arguments to pass to the callable

        :param run_after: the futures representing other callables that must complete successfully before
                          this callable is invoked

        :return: a Future instance representing the callable
        """
        return DeferredFuture(self, key, run_after, callable_, *args, **kwargs)

    def _submit(self, deferred_future: DeferredFuture) -> Future:
        future = self._tpe.submit(*deferred_future.args, **deferred_future.kwargs)
        self._futures[future] = deferred_future.key
        return future

    def _collect_futures(self) -> ErrorsByKey:
        errors = defaultdict(list)
        while self._futures:
            for future in as_completed(self._futures.keys()):
                key = self._futures.pop(future)
                e = future.exception()
                if e is not None:
                    errors[key].append(e)
        return errors

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._errors = self._collect_futures()
        self._tpe.shutdown(wait=True)
        return False

    @property
    def errors(self):
        assert self._errors is not None, 'Errors can only be used after context manager has been exited.'
        return self._errors


def test():
    expected = 0
    actual = 0
    with DeferredTaskExecutor() as executor:
        def g(n, i):
            nonlocal actual, expected
            actual += 1
            print(n, i)
            time.sleep(random() * 10)
            if n:
                m = int(1 << n)
                expected += m
                parents = [executor.defer(m + i, g, n - 1, j) for j in range(m)]
                if False and i and random() < .5:
                    raise RuntimeError(n, i)
                expected += m
                [executor.defer(m + i, g, n - 1, -j, run_after=parents) for j in range(m)]

        g(3, 0)
    print(executor.errors)
    print(expected, actual)


if __name__ == '__main__':
    test()
