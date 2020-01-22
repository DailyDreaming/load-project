import os
import pathlib


# TODO: https://codereview.stackexchange.com/questions/162426/subclassing-pathlib-path

class Path(pathlib.PosixPath):

    # _flavour = pathlib.Path._flavour

    # Work around https://bugs.python.org/issue30618, fixed on 3.7+

    def readlink(self) -> 'Path':
        """
        Return the path to which the symbolic link points.
        """
        path = self._accessor.readlink(self)
        obj = self._from_parts((path,), init=False)
        obj._init(template=self)
        return obj

    def follow(self) -> 'Path':
        """
        This methods performs one level of symbolic link resolution. For paths
        representing a symbolic link with an absolute target, this methods is
        equivalent to readlink(). For symbolic links with relative targets, this
        method returns the result of appending the target to the parent of this
        path. Unless you need the target of the symbolic link verbatim, you
        should prefer this method over readlink().
        """
        target = self.readlink()
        if target.is_absolute():
            return target
        else:
            return self.parent / target

    # Sorely needed, added in 3.8

    def link_to(self, target: 'Path'):
        """
        Create a hard link pointing to a path named target.
        """
        if self._closed:
            self._raise_closed()
        os.link(str(self), str(target))
