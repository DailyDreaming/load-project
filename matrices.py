import logging
import shutil
import sys
from abc import abstractmethod, ABCMeta
from pathlib import Path

import os
from typing import NamedTuple

log = logging.getLogger(__file__)


class Matrix(NamedTuple):
    mtx: str
    genes: str
    barcodes: str


# noinspection PyUnusedLocal
def convert_csv_to_mtx(input_file: Path, output_dir: Path, delimiter: str = ',', rows_are_genes: bool = True):
    raise NotImplementedError()


class Converter(metaclass=ABCMeta):

    def __init__(self, project_dir: Path):
        self.project_dir = project_dir

    @property
    def matrices_dir(self) -> Path:
        return self.project_dir / 'matrices'

    @property
    def geo_dir(self) -> Path:
        return self.project_dir / 'geo'

    @property
    def bundle_dir(self) -> Path:
        return self.project_dir / 'bundle'

    @property
    def zip_file(self) -> Path:
        return self.bundle_dir / 'matrix.mtx.zip'

    def matrix_dir(self, input_: str) -> Path:
        return self.matrices_dir / input_.replace('/', '__')

    def convert(self):
        if self.zip_file.exists():
            log.info('Final matrix already exists for project %s; moving on.', self.project_dir)
        else:
            self._convert()
            self._create_zip()

    def _create_zip(self):
        final_matrix = self.zip_file
        os.makedirs(final_matrix.parent, exist_ok=True)
        # make_archive adds it's own .zip at the end
        shutil.make_archive(self.zip_file.parent / self.zip_file.stem, 'zip', self.matrices_dir)

    @abstractmethod
    def _convert(self):
        raise NotImplementedError()

    def _link_matrix(self, matrix: 'Matrix'):
        matrix_dir = self.matrix_dir(matrix.mtx)
        matrix_dir.mkdir(parents=True, exist_ok=True)
        assert all(name.endswith('.gz') for name in matrix)
        idempotent_link(self.geo_dir / matrix.mtx, matrix_dir / 'matrix.mtx.gz')
        idempotent_link(self.geo_dir / matrix.genes, matrix_dir / 'genes.tsv.gz')
        idempotent_link(self.geo_dir / matrix.barcodes, matrix_dir / 'barcodes.tsv.gz')

    def _link_matrices(self, matrices):
        for matrix in matrices:
            self._link_matrix(matrix)


def idempotent_link(src: Path, dst: Path):
    if src.stat().st_ino != dst.stat().st_ino:
        dst.unlink()
        os.link(src, dst)


class GSE107909(Converter):
    """
    04ba7269-1301-5758-8f13-025565326f66
    """

    def _convert(self):
        csvs = [
            "GSE107909_RAW/GSM2883183_PLNr9c.csv.gz",
            "GSE107909_RAW/GSM2883182_PLN++.csv.gz"
        ]
        for csv in csvs:
            convert_csv_to_mtx(self.geo_dir / csv, self.matrix_dir(csv), ',', rows_are_genes=True)


class GSE106273(Converter):
    """
    099c02da-23b2-5748-8618-92bc6770dc51
    """

    def _convert(self):
        self._link_matrices([
            Matrix(
                mtx='GSE106273_combined_matrix.tsv.gz',
                genes='GSE106273_combined_genes.tsv.gz',
                barcodes='GSE106273_combined_barcodes.tsv.gz',
            ),
        ])


class GSE117089(Converter):
    """
    061ec9d5-9acf-54db-9eee-555136d5ce41
    """

    def _convert(self):
        self._link_matrices([
            Matrix(
                barcodes='GSE117089_RAW/GSM3271042_RNA_only_A549_cell.txt.gz',
                genes='GSE117089_RAW/GSM3271042_RNA_only_A549_gene.txt.gz',
                mtx='GSE117089_RAW/GSM3271042_RNA_only_A549_gene_count.txt.gz',
            ),
            Matrix(
                barcodes='GSE117089_RAW/GSM3271043_ATAC_only_A549_cell.txt.gz',
                genes='GSE117089_RAW/GSM3271043_ATAC_only_A549_peak.txt.gz',
                mtx='GSE117089_RAW/GSM3271043_ATAC_only_A549_peak_count.txt.gz',
            ),
            Matrix(
                barcodes='GSE117089_RAW/GSM3271044_RNA_mouse_kidney_cell.txt.gz',
                genes='GSE117089_RAW/GSM3271044_RNA_mouse_kidney_gene.txt.gz',
                mtx='GSE117089_RAW/GSM3271044_RNA_mouse_kidney_gene_count.txt.gz',
            ),
            Matrix(
                barcodes='GSE117089_RAW/GSM3271045_ATAC_mouse_kidney_cell.txt.gz',
                genes='GSE117089_RAW/GSM3271045_ATAC_mouse_kidney_peak.txt.gz',
                mtx='GSE117089_RAW/GSM3271045_ATAC_mouse_kidney_peak_count.txt.gz',
            )
        ])


def synthesize_matrices(projects: Path):
    failed_projects = {}
    succeeded_projects = set()
    for project_dir in projects.iterdir():
        if project_dir.is_symlink():
            try:
                converter_class = globals()[project_dir.name]
                converter = converter_class(project_dir)
                converter.convert()
            except Exception as e:
                failed_projects[project_dir] = e
                log.exception('Failed to process project', exc_info=True)
            else:
                succeeded_projects.add(project_dir)

    print('Failed projects', file=sys.stderr)
    for p in failed_projects:
        print(p, file=sys.stderr)

    print('Succeeded projects')
    for project_dir in succeeded_projects:
        print(project_dir)


if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s',
                        level=logging.DEBUG)

    # To download some projects for testing run:
    # scp -r ubuntu@skunk.dev.explore.data.humancellatlas.org:/home/ubuntu/load-project/projects/0* ./test/projects
    input_dir = sys.argv[1]
    synthesize_matrices(Path(input_dir))
