import os
from pathlib import Path
from typing import (
    Tuple,
    Sequence,
)

from convert_csv import geo_dir


def one_mtx(project_dir: Path, filenames: Tuple[str, str, str]):
    """
    Only a single mtx triplet with incorrect names
    """
    geo = geo_dir(project_dir)
    for src, dst in zip(filenames, ('genes.tsv.gz', 'barcodes.tsv.gz', 'matrix.mtx.gz')):
        try:
            os.link(str(geo / src), str(project_dir / dst))
        except FileExistsError:
            pass

