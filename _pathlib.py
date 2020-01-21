import os
import pathlib


# TODO: https://codereview.stackexchange.com/questions/162426/subclassing-pathlib-path

class Path(pathlib.PosixPath):

    # _flavour = pathlib.Path._flavour

    # Work around https://bugs.python.org/issue30618, fixed on 3.7+

    def readlink(self):
        """
        Return the path to which the symbolic link points.
        """
        path = self._accessor.readlink(self)
        obj = self._from_parts((path,), init=False)
        obj._init(template=self)
        return obj

    # Sorely needed, added in 3.8

    def link_to(self, target: 'Path'):
        """
        Create a hard link pointing to a path named target.
        """
        if self._closed:
            self._raise_closed()
        os.link(str(self), str(target))
