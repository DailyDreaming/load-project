"""
Convert a 10xgenomics HDF5 matrix to MTX.

https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices
http://scipy-lectures.org/advanced/scipy_sparse/csc_matrix.html
https://math.nist.gov/MatrixMarket/formats.html
"""
from itertools import chain
import logging
from pathlib import Path
from typing import (
    Any,
    Iterable,
    Tuple,
)

from dataclasses import (
    dataclass,
    fields,
)
import h5py
from h5py import (
    Dataset,
    Group,
)
from more_itertools import (
    one,
    pairwise,
)

from csv2mtx import write_gzip_file

log = logging.getLogger(__name__)


@dataclass
class Matrix:
    file_path: Path
    name: str
    barcodes: Dataset
    data: Dataset
    gene_names: Dataset
    genes: Dataset
    indices: Dataset
    indptr: Dataset
    shape: Dataset

    @classmethod
    def from_group(cls, file_path: Path, group: Group) -> 'Matrix':
        return cls(file_path, group.name, *(group[f.name] for f in fields(cls)[2:]))

    @property
    def rows(self) -> int:
        return self.shape[0]

    @property
    def columns(self) -> int:
        return self.shape[1]

    @property
    def cells(self) -> int:
        return len(self.data)

    def to_mtx(self, output_dir: Path):
        assert self.rows == len(self.genes)
        assert self.columns == len(self.barcodes)
        assert len(self.data) == len(self.indices)
        write_gzip_file(output_dir / 'genes.gz', chain(['genes'], map(bytes.decode, self.genes)))
        write_gzip_file(output_dir / 'barcodes.gz', chain(['barcodes'], map(bytes.decode, self.barcodes)))
        log.info('Reading matrix data from %s ...', self.file_path)
        write_gzip_file(output_dir / 'matrix.mtx.gz', chain([
            '%%MatrixMarket matrix coordinate integer general',
            f'{self.rows} {self.columns} {self.cells}'
        ], self._mtx_lines()))

    def _mtx_lines(self) -> Iterable[str]:
        # Sorting shouldn't be necessary since the MatrixMarket format does not
        # dictate an ordering of lines in the .mtx but most of the .mtx files
        # I've seen are sorted by row index first and column index second. If
        # we didn't sort, we'd get the transposed ordering by column first and
        # row second. This is because the HDF5 matrices are column-oriented.
        return (' '.join(map(str, t)) for t in sorted(self._mtx_tuples()))

    def _mtx_tuples(self) -> Iterable[Tuple[int, int, Any]]:
        for column, (start, end) in enumerate(pairwise(self.indptr)):
            data = self.data[start:end]
            rows = self.indices[start:end]
            for row, value in zip(rows, data):
                yield row + 1, column + 1, value


def convert_h5_to_mtx(input_file: Path, output_dir: Path) -> None:
    with h5py.File(str(input_file), mode='r') as h5:
        group = one(h5.values())
        m = Matrix.from_group(input_file, group)
        output_dir.mkdir(parents=True, exist_ok=True)  # FIXME: move to convert_matrices.py
        m.to_mtx(output_dir)
