import argparse
import csv
import json
import logging
import sys
from typing import (
    MutableMapping,
    Sequence,
    Tuple,
    Union,
)

from _pathlib import Path
from util import (
    get_target_project_dirs,
    open_maybe_gz,
    update_project_stats,
)


class CountCells:

    def __init__(self, argv):
        logging.basicConfig(level=logging.INFO)
        parser = argparse.ArgumentParser(description=__doc__)
        group = parser.add_mutually_exclusive_group(required=True)
        group.add_argument('--list', '-l',
                           help='List cell count values',
                           action='store_true')
        group.add_argument('--write', '-w',
                           help='Write cell count files',
                           action='store_true')
        parser.add_argument('--verbose', '-v',
                            action='store_true',
                            help='Verbose debug output')
        parser.add_argument('--slow', '-s',
                            action='store_true',
                            help='Get cell count by counting unique barcode/gene indexes.')
        parser.add_argument('--medium', '-m',
                            action='store_true',
                            help='Get cell count by counting lines in the barcode/gene file.')
        parser.add_argument('--fast', '-f',
                            action='store_true',
                            help='Get cell count from the header of the matrix file.')
        self.args = parser.parse_args(argv)

        if not (self.args.slow or self.args.medium or self.args.fast):
            parser.error('Error: One or more [--slow, --medium, --fast] argument must be specified.')

    def run(self):
        if self.args.verbose:
            logging.getLogger().setLevel(logging.DEBUG)

        count_method = []
        if self.args.slow:
            count_method.append('slow')
        if self.args.medium:
            count_method.append('medium')
        if self.args.fast:
            count_method.append('fast')

        for project_dir in get_target_project_dirs():
            cell_count, gene_count = self.get_project_contents_count(project_dir, count_method)
            if self.args.write:
                self.write_cell_gene_count(project_dir, cell_count, gene_count)

    @classmethod
    def get_cached_cell_count(cls, project_dir: Path) -> Union[int, None]:
        """
        Return a cell count for the project with the given accession.
        """
        stats_file = project_dir / 'stats.json'
        if stats_file.exists():
            with open(str(stats_file), 'r') as f:
                cell_count = json.load(f).get('cell_count', None)
            return cell_count
        else:
            return 0

    @classmethod
    def get_cached_gene_count(cls, project_dir: Path) -> Union[int, None]:
        """
        Return a gene count for the project with the given accession.
        """
        stats_file = project_dir / 'stats.json'
        if stats_file.exists():
            with open(str(stats_file), 'r') as f:
                gene_count = json.load(f).get('gene_count', None)
            return gene_count
        else:
            return 0

    @classmethod
    def get_cached_cell_counts(cls) -> MutableMapping[str, int]:
        """
        Return a mapping from accessions to cell counts.
        """

        return {
            p.name: cls.get_cached_cell_count(p)
            for p in get_target_project_dirs()
        }

    @classmethod
    def get_cached_gene_counts(cls) -> MutableMapping[str, int]:
        """
        Return a mapping from accessions to gene counts.
        """

        return {
            p.name: cls.get_cached_gene_count(p)
            for p in get_target_project_dirs()
        }

    @classmethod
    def write_cell_gene_count(cls, project_dir: Path, cell_count: int, gene_count: int):
        """
        Write the accession cell and gene counts to the project's stats JSON file.
        """

        with update_project_stats(project_dir) as stats:
            stats['cell_count'] = cell_count
            stats['gene_count'] = gene_count

    def get_project_contents_count(self, project_dir: Path, count_method: Sequence[str]) -> Tuple[int, int]:
        """
        Count the number of cells in a project.

        :param project_dir: Path to the project directory
        :param count_method: One or more counting methods to use ('slow', 'medium', 'fast')
        """

        cell_counts = {key: 0 for key in count_method}
        gene_counts = {key: 0 for key in count_method}

        matrix_dir = project_dir / 'matrices'
        for mtx_file in matrix_dir.glob('**/matrix.mtx.gz'):
            cell_count = {}
            gene_count = {}
            logging.debug('Getting counts from %s', mtx_file.parent)
            if 'slow' in count_method:
                cell_count['slow'] = self.count_unique_index_values(mtx_file, col=2)
                cell_counts['slow'] += cell_count['slow']
                gene_count['slow'] = self.count_unique_index_values(mtx_file, col=1)
                gene_counts['slow'] += gene_count['slow']
            if 'medium' in count_method:
                barcode_file = mtx_file.parent / 'barcodes.tsv.gz'
                cell_count['medium'] = self.count_rows_in_tsv(barcode_file)
                cell_counts['medium'] += cell_count['medium']
                gene_file = mtx_file.parent / 'genes.tsv.gz'
                gene_count['medium'] = self.count_rows_in_tsv(gene_file)
                gene_counts['medium'] += gene_count['medium']
            if 'fast' in count_method:
                cell_count['fast'] = self.get_count_from_mtx_header(mtx_file, col=2)
                cell_counts['fast'] += cell_count['fast']
                gene_count['fast'] = self.get_count_from_mtx_header(mtx_file, col=1)
                gene_counts['fast'] += gene_count['fast']

            for speed in ('slow', 'medium', 'fast'):
                if speed in count_method:
                    logging.info('%s cell count of %s is %s',
                                 speed.capitalize(), barcode_file if speed == 'medium' else mtx_file, cell_count[speed])
            for speed in ('slow', 'medium', 'fast'):
                if speed in count_method:
                    logging.info('%s gene count of %s is %s',
                                 speed.capitalize(), gene_file if speed == 'medium' else mtx_file, gene_count[speed])
            if len(cell_count) > 1 and len(set(cell_count.values())) > 1:
                logging.warning('Cell count mismatch in %s: %s', mtx_file.parent, cell_count)
            if len(gene_count) > 1 and len(set(gene_count.values())) > 1:
                logging.warning('Gene count mismatch in %s: %s', mtx_file.parent, gene_count)

        total_cell_count = max(cell_counts.values())  # TODO: prioritize value from one method over another?
        total_gene_count = max(gene_counts.values())
        logging.info('Total cell count in %s is %s', project_dir, total_cell_count)
        logging.info('Total gene count in %s is %s', project_dir, total_gene_count)

        # Compare counts between counting methods
        if len(count_method) > 1:
            if len(set(cell_counts.values())) > 1:
                logging.warning('Cell count mismatch in %s: %s', project_dir, cell_counts)
            if len(set(gene_counts.values())) > 1:
                logging.warning('Gene count mismatch in %s: %s', project_dir, gene_counts)

        return total_cell_count, total_gene_count

    @classmethod
    def sniff_for_dialect(cls, matrix_file: Path) -> csv.Sniffer:
        with open_maybe_gz(matrix_file, 'rt', newline='') as csv_file:
            first_line_length = len(csv_file.readline().strip())
            assert first_line_length > 0, f'File has no first line "{matrix_file}"'
            csv_file.seek(0)
            return csv.Sniffer().sniff(csv_file.read(first_line_length))

    def count_unique_index_values(self, matrix_file: Path, col: int) -> int:
        indexes = None
        dialect = self.sniff_for_dialect(matrix_file)
        with open_maybe_gz(matrix_file, 'rt', newline='') as csv_file:
            csv_reader = csv.reader(csv_file, dialect)
            for line_num, row in enumerate(csv_reader):
                if not row[0].startswith('%'):  # skip comment lines in the csv
                    assert len(row) == 3, f'{matrix_file} has line {line_num} with {len(row)} columns instead of 3'
                    if indexes is None:  # skip first non-comment line (a header line)
                        indexes = set()
                    else:
                        indexes.add(row[col - 1])
        return len(indexes)

    def get_count_from_mtx_header(self, matrix_file: Path, col: int) -> int:
        dialect = self.sniff_for_dialect(matrix_file)
        with open_maybe_gz(matrix_file, 'rt', newline='') as csv_file:
            csv_reader = csv.reader(csv_file, dialect)
            for line_num, row in enumerate(csv_reader):
                if not row[0].startswith('%'):  # skip comment lines in the csv
                    assert len(row) == 3, f'{matrix_file} has line {line_num} with {len(row)} columns instead of 3'
                    # The first non-comment line contains the dimensions of the
                    # matrix and the number of non-zero cells in that matrix.
                    return int(row[col - 1])

    def count_rows_in_tsv(self, file: Path) -> int:
        # We're not using the `csv` module to save time, avoiding to split each
        # row and create a list per row.
        with open_maybe_gz(file, 'rt') as f:
            header_row = next(f)
            num_cols = header_row.count('\t') + 1
            rows = 0
            for row in f:
                num_cells = row.count('\t') + 1
                rows += 1  # increment before assert so that error message report the correct line number
                assert num_cols == num_cells, f'{file} has line {rows} with {num_cells} column(s) instead of {num_cols}'
            return rows


if __name__ == '__main__':
    CountCells(sys.argv[1:]).run()
