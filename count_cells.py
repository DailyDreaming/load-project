from collections import Sequence

import argparse
import csv
import json
import logging
import os
from pathlib import Path
import sys
from typing import (
    Union,
    Sequence,
)

from util import (
    open_maybe_gz,
    get_target_project_dirs,
)

cell_counts_file = Path('cell_counts.json')


def main(argv):
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--list', '-l',
                       help='list cell count values',
                       action='store_true')
    group.add_argument('--write', '-w',
                       help='write cell count files',
                       action='store_true')
    parser.add_argument('--verbose', '-v',
                        action='store_true',
                        help='Verbose debug output')
    args = parser.parse_args(argv)

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    project_dirs = get_target_project_dirs()

    for project_dir in project_dirs:
        count_one_project_cells(project_dir.name, args.write)


def count_one_project_cells(accession_id: str, save_to_file=True):
    """
    Count cells in one project

    :param accession_id: An accession id with a downloaded matrix file
    :param save_to_file: If true write the cell count total(s) to a global cell count file
    """
    cell_count = get_project_cell_count(accession_id)
    if save_to_file:
        update_cell_count_file(accession_id, cell_count)


def get_accession_ids() -> Sequence[str]:
    """
    Return a list of the accession ids
    """
    accession_ids = []
    for p in Path('projects').iterdir():
        logging.debug('Checking: %s', p)
        if p.is_dir() and p.is_symlink():
            accession_ids.append(p.name)
            logging.debug('Found: %s', p.name)
    return accession_ids


def get_cell_counts():
    if cell_counts_file.exists():
        with open(cell_counts_file, 'r') as f:
            cell_counts = json.loads(f.read())
        return cell_counts
    return {}


def update_cell_count_file(accession_id: str, cell_count: int) -> bool:
    """
    Write the accession cell count to the global cell count file

    :param accession_id: An accession id
    :param cell_count: An int value to save, or None to remove entry
    :return: Boolean status of the save
    """
    if not accession_id:
        return False
    cell_counts = get_cell_counts()
    if isinstance(cell_count, int):
        logging.info('Writing accession %s cell count.', accession_id)
        cell_counts[accession_id] = cell_count
    elif cell_count is None and accession_id in cell_counts:
        logging.info('Removing accession %s cell count.', accession_id)
        del cell_counts[accession_id]
    with open(cell_counts_file, 'w') as f:
        f.write(json.dumps(cell_counts, sort_keys=True, indent='    '))
    return True


def get_project_cell_count(accession_id: str) -> Union[int, None]:
    """
    Get the cell count from a project's matrix file(s)

    :param accession_id: An accession id that has downloaded matrix file(s)
    :return: A count of cells
    """
    if not os.path.isdir(f'projects/{accession_id}'):
        logging.warning('Unable to find path %s', f'projects/{accession_id}')
        return None
    total_count = 0
    for file in Path(f'projects/{accession_id}/matrices').glob('**/*'):
        if file.is_file() and file.name.endswith(('.mtx', '.mtx.gz')):
            cell_count = count_cells(file)
            logging.info('Cell count in %s is %s', file, cell_count)
            total_count += cell_count
    logging.info('Total cell count in %s is %s', accession_id, total_count)
    return total_count


def count_cells(matrix_file: Path) -> Union[int, None]:
    """
    Count the number of cells in a matrix file

    :param matrix_file: A mtx file with tab/space/comma delimited data
    :return: A count of cells found in the given file
    """
    if matrix_file.exists():
        with open_maybe_gz(matrix_file, 'rt', newline='') as csv_file:
            first_line_length = len(csv_file.readline().strip())
            csv_file.seek(0)
            if first_line_length:
                dialect = csv.Sniffer().sniff(csv_file.read(first_line_length))
            else:
                logging.error(f'File has no first line "{matrix_file}"')
                return None
        barcode_indexes = set()
        with open_maybe_gz(matrix_file, 'rt', newline='') as csv_file:
            csv_reader = csv.reader(csv_file, dialect)
            for row in csv_reader:
                if row[0].startswith('%'):
                    continue  # skip comment lines in the csv
                if len(row) == 3:
                    barcode_indexes.add(row[1])  # 2nd column is barcode
        return len(barcode_indexes)
    else:
        logging.error(f'File not found "{matrix_file}"')
        return None


if __name__ == '__main__':
    main(sys.argv[1:])
