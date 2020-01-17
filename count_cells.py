import argparse
import csv
import json
import logging
from _pathlib import Path
import sys
from typing import (
    Sequence,
    Mapping,
    Union,
)

from util import (
    generate_project_uuid,
    get_target_project_dirs,
    open_maybe_gz,
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
        parser.add_argument('--accession', '-a', type=str,
                            help='Count cells for a single accession '
                                 'instead of all')
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

        if self.args.accession:
            accessions = [self.args.accession]
        else:
            accessions = [p.name for p in get_target_project_dirs()]

        for accession in accessions:
            cell_count = self.get_project_cell_count(accession)
            if self.args.write:
                self.write_cell_count(accession, cell_count)

    @classmethod
    def get_cached_cell_count(cls, accession: str) -> int:
        """
        Return a cell count for the given accession
        """
        stats_file = Path(f'projects/{accession}/stats.json')
        if stats_file.exists():
            with open(stats_file, 'r') as f:
                cell_count = json.loads(f.read()).get('cell_count', 0)
            return cell_count
        else:
            return 0

    @classmethod
    def get_cached_cell_counts(cls) -> Mapping[str, int]:
        """
        Return a dict of accession id to cell count
        """
        cell_counts = {}
        for accession in [p.name for p in get_target_project_dirs()]:
            cell_counts[accession] = cls.get_cached_cell_count(accession)
        return cell_counts

    @classmethod
    def write_cell_count(cls, accession: str, cell_count: int) -> bool:
        """
        Write the accession cell count to the project's stats json

        :param accession: An accession id
        :param cell_count: An int value to write to the json file
        :return: Status of the write
        """
        stats_file = Path(f'projects/{accession}/stats.json')
        if not stats_file.exists():
            with open(stats_file, 'w') as f:
                logging.info('Writing %s', stats_file)
                stats = {
                    'project_uuid': generate_project_uuid(accession),
                    'cell_count': cell_count
                }
                f.write(json.dumps(stats, sort_keys=True, indent='    '))
        else:
            with open(stats_file, 'r+') as f:
                logging.info('Updating %s', stats_file)
                stats = json.loads(f.read())
                stats['cell_count'] = cell_count
                f.seek(0)
                f.write(json.dumps(stats, sort_keys=True, indent='    '))
        return True

    def get_project_cell_count(self, accession: str) -> Union[int, None]:
        """
        Count the number of cells in a project

        :param accession: An accession id that has downloaded matrix file(s)
        :return: A count of cells
        """
        total_cell_count = 0
        for mtx_file in Path(f'projects/{accession}/matrices').glob('**/matrix.mtx.gz'):
            cell_count_from_matrix = self.count_cells(mtx_file)
            logging.info('Cell count in %s is %s', mtx_file, cell_count_from_matrix)
            if cell_count_from_matrix is not None:
                total_cell_count += cell_count_from_matrix
            barcodes_file = mtx_file.parent / 'barcodes.tsv.gz'
            cell_count_from_barcodes = self.count_cells_from_barcodes(barcodes_file)
            if cell_count_from_matrix != cell_count_from_barcodes:
                logging.warning('Cell count mismatch found for %s: %s vs %s',
                                mtx_file.parent, cell_count_from_matrix, cell_count_from_barcodes)
        logging.info('Total cell count in %s is %s', accession, total_cell_count)
        return total_cell_count

    def count_cells(self, matrix_file: Path) -> Union[int, None]:
        """
        Count the number of cells in a matrix file

        :param matrix_file: A mtx file with tab/space/comma delimited data
        :return: A count of cells found in the given file
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
