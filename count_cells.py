import argparse
import csv
import logging
import json
from _pathlib import Path
import sys
from typing import (
    MutableMapping,
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
        parser.add_argument('--fast', '-f',
                            action='store_true',
                            help='Use barcodes dimension from MTX header '
                                 'instead of counting unique barcode indices.')
        self.args = parser.parse_args(argv)

    def run(self):
        if self.args.verbose:
            logging.getLogger().setLevel(logging.DEBUG)

        for project_dir in get_target_project_dirs():
            cell_count = self.get_project_cell_count(project_dir)
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

    def get_project_cell_count(self, project_dir: Path) -> int:
        """
        Count the number of cells in a project.
        """

        total_cell_count = 0
        matrix_dir = project_dir / 'matrices'
        for mtx_file in matrix_dir.glob('**/matrix.mtx.gz'):
            cell_count_from_matrix = self.count_cells(mtx_file)
            logging.info('Cell count in %s is %s', mtx_file, cell_count_from_matrix)
            if cell_count_from_matrix is not None:
                total_cell_count += cell_count_from_matrix
            barcodes_file = mtx_file.parent / 'barcodes.tsv.gz'
            cell_count_from_barcodes = self.count_cells_from_barcodes(barcodes_file)
            if cell_count_from_matrix != cell_count_from_barcodes:
                logging.warning('Cell count mismatch found for %s: %s vs %s',
                                mtx_file.parent, cell_count_from_matrix, cell_count_from_barcodes)
        logging.info('Total cell count in %s is %s', project_dir, total_cell_count)
        return total_cell_count

    def count_cells(self, matrix_file: Path) -> Union[int, None]:
        """
        Count the number of cells in a matrix file.

        :param matrix_file: A mtx file with tab/space/comma delimited data.
        :return: A count of cells found in the given file.
        """

        with open_maybe_gz(matrix_file, 'rt', newline='') as csv_file:
            first_line_length = len(csv_file.readline().strip())
            assert first_line_length > 0, f'File has no first line "{matrix_file}"'
            csv_file.seek(0)
            dialect = csv.Sniffer().sniff(csv_file.read(first_line_length))
        barcode_indexes = None
        with open_maybe_gz(matrix_file, 'rt', newline='') as csv_file:
            csv_reader = csv.reader(csv_file, dialect)
            for line_num, row in enumerate(csv_reader):
                if not row[0].startswith('%'):  # skip comment lines in the csv
                    assert len(row) == 3, f'{matrix_file} has line {line_num} with {len(row)} columns instead of 3'
                    if barcode_indexes is None:
                        # The first non-comment line contains the dimensions
                        # of the matrix and the number of non-zero cells in
                        # that matrix.
                        if self.args.fast:
                            return int(row[1])
                        else:
                            barcode_indexes = set()
                    else:
                        barcode_indexes.add(row[1])  # 2nd column is barcode
        return len(barcode_indexes)

    def count_cells_from_barcodes(self, barcodes_file: Path) -> Union[int, None]:
        with open_maybe_gz(barcodes_file,  'rt') as f:
            next(f)  # skip header line
            return sum(1 for _ in f)


if __name__ == '__main__':
    CountCells(sys.argv[1:]).run()
