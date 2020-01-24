import argparse
import csv
import logging
import json
from _pathlib import Path
import sys
from typing import (
    MutableMapping,
    Sequence,
    Union,
)
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
                            help='Get cell count by counting unique barcode indexes.')
        parser.add_argument('--medium', '-m',
                            action='store_true',
                            help='Get cell count by counting lines in the barcodes file.')
        parser.add_argument('--fast', '-f',
                            action='store_true',
                            help='Get cell count from the barcodes dimension in the MTX header.')
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
            cell_count = self.get_project_cell_count(project_dir, count_method)
            if self.args.write:
                self.write_cell_count(project_dir, cell_count)

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
    def get_cached_cell_counts(cls) -> MutableMapping[str, int]:
        """
        Return a mapping from accessions to cell counts.
        """

        return {
            p.name: cls.get_cached_cell_count(p)
            for p in get_target_project_dirs()
        }

    @classmethod
    def write_cell_count(cls, project_dir: Path, cell_count: int):
        """
        Write the accession cell count to the project's stats JSON file.
        """

        with update_project_stats(project_dir) as stats:
            stats['cell_count'] = cell_count

    def get_project_cell_count(self, project_dir: Path, count_method: Sequence[str]) -> int:
        """
        Count the number of cells in a project.

        :param project_dir: Path to the project directory
        :param count_method: One or more counting methods to use ('slow', 'medium', 'fast')
        """

        cell_counts = {key: 0 for key in count_method}

        matrix_dir = project_dir / 'matrices'
        for mtx_file in matrix_dir.glob('**/matrix.mtx.gz'):
            cell_count = {}
            if 'slow' in count_method:
                cell_count['slow'] = self.count_cells_in_matrix(mtx_file)
                logging.info('Slow cell count of %s is %s', mtx_file, cell_count['slow'])
                cell_counts['slow'] += cell_count['slow']
            if 'medium' in count_method:
                barcodes_file = mtx_file.parent / 'barcodes.tsv.gz'
                cell_count['medium'] = self.count_cells_from_barcodes(barcodes_file)
                logging.info('Medium cell count of %s is %s', barcodes_file, cell_count['medium'])
                cell_counts['medium'] += cell_count['medium']
            if 'fast' in count_method:
                cell_count['fast'] = self.count_cells_in_mtx_header(mtx_file)
                logging.info('Fast cell count of %s is %s', mtx_file, cell_count['fast'])
                cell_counts['fast'] += cell_count['fast']
            if len(cell_count) > 1 and len(set(cell_count.values())) > 1:
                logging.warning('Count mismatch in %s: %s',
                                mtx_file.parent, cell_count)

        total_count = max(cell_counts.values())  # TODO: prioritize value from one method over another?
        logging.info('Total cell count in %s is %s', project_dir, total_count)

        # Compare counts between counting methods
        if len(count_method) > 1 and len(set(cell_counts.values())) > 1:
            logging.warning('Count mismatch in %s: %s', project_dir, cell_counts)

        return total_count

    @classmethod
    def sniff_for_dialect(cls, matrix_file: Path) -> csv.Sniffer:
        with open_maybe_gz(matrix_file, 'rt', newline='') as csv_file:
            first_line_length = len(csv_file.readline().strip())
            assert first_line_length > 0, f'File has no first line "{matrix_file}"'
            csv_file.seek(0)
            return csv.Sniffer().sniff(csv_file.read(first_line_length))

    def count_cells_in_matrix(self, matrix_file: Path) -> int:
        barcode_indexes = None
        dialect = self.sniff_for_dialect(matrix_file)
        with open_maybe_gz(matrix_file, 'rt', newline='') as csv_file:
            csv_reader = csv.reader(csv_file, dialect)
            for line_num, row in enumerate(csv_reader):
                if not row[0].startswith('%'):  # skip comment lines in the csv
                    assert len(row) == 3, f'{matrix_file} has line {line_num} with {len(row)} columns instead of 3'
                    if barcode_indexes is None:  # skip first non-comment line (a header line)
                        barcode_indexes = set()
                    else:
                        barcode_indexes.add(row[1])  # 2nd column is barcode
        return len(barcode_indexes)

    def count_cells_in_mtx_header(self, matrix_file: Path) -> int:
        dialect = self.sniff_for_dialect(matrix_file)
        with open_maybe_gz(matrix_file, 'rt', newline='') as csv_file:
            csv_reader = csv.reader(csv_file, dialect)
            for line_num, row in enumerate(csv_reader):
                if not row[0].startswith('%'):  # skip comment lines in the csv
                    assert len(row) == 3, f'{matrix_file} has line {line_num} with {len(row)} columns instead of 3'
                    # The first non-comment line contains the dimensions of the
                    # matrix and the number of non-zero cells in that matrix.
                    return int(row[1])

    def count_cells_from_barcodes(self, barcodes_file: Path) -> int:
        with open_maybe_gz(barcodes_file,  'rt') as f:
            next(f)  # skip header line
            return sum(1 for _ in f)


if __name__ == '__main__':
    CountCells(sys.argv[1:]).run()
