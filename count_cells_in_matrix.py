import argparse
import csv
import json
import os
import sys

from typing import (
    Sequence,
    Union
)

import logging
logging.basicConfig(level=logging.INFO)


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--list', '-l',
                       metavar='ACCESSION',
                       help='list cell count value for given accession id')
    group.add_argument('--write', '-w',
                       metavar='ACCESSION',
                       help='write a cell count file for given accession id')
    group.add_argument('--list-all', '-L',
                       action='store_true',
                       help='list cell counts for all accession ids')
    group.add_argument('--write-all', '-W',
                       action='store_true',
                       help='write cell count files for all accession ids')
    parser.add_argument('--verbose', '-v',
                        action='store_true',
                        help='Verbose debug output')
    args = parser.parse_args(argv)

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    if args.list:
        count_one_project_cells(args.list, save_to_file=False)
    elif args.write:
        count_one_project_cells(args.write, save_to_file=True)
    elif args.list_all:
        count_all_project_cells(save_to_file=False)
    elif args.write_all:
        count_all_project_cells(save_to_file=True)


def count_all_project_cells(save_to_file=True):
    """
    Count cells in all projects

    :param save_to_file: If true write the cell count totals to a global cell count file
    """
    for accession_id in get_accession_ids():
        logging.info('Checking: %s ...', accession_id)
        count_one_project_cells(accession_id, save_to_file)


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
    for item in os.listdir('projects'):
        path = f'projects/{item}'
        logging.debug('Checking: %s', path)
        if os.path.isdir(path) and os.path.islink(path):
            accession_ids.append(item)
    return accession_ids


def update_cell_count_file(accession_id: str, cell_count: int) -> bool:
    """
    Write the accession cell count to the global cell count file

    :param accession_id: An accession id
    :param cell_count: An int value to save, or None to remove entry
    :return: Boolean status of the save
    """
    if not accession_id:
        return False
    cell_counts_file = 'cell_counts.json'
    if os.path.exists(cell_counts_file):
        with open(cell_counts_file, 'r') as f:
            cell_counts = json.loads(f.read())
    else:
        cell_counts = {}
    if isinstance(cell_count, int):
        logging.info('Writing accession %s cell count.', accession_id)
        cell_counts[accession_id] = cell_count
    elif cell_count is None and accession_id in cell_counts:
        logging.info('Removing accession %s cell count.', accession_id)
        del(cell_counts[accession_id])
    with open(cell_counts_file, 'w') as f:
        f.write(json.dumps(cell_counts, sort_keys=True, indent='    '))
    return True


def get_project_cell_count(accession_id: str) -> Union[int, None]:
    """
    Get the cell count from a project's matrix file

    :param accession_id: An accession id that has a downloaded matrix file
    :return: A count of cells
    """
    matrix_file_full = f'projects/{accession_id}/matrix.mtx'
    if not os.path.isdir(f'projects/{accession_id}') or not os.path.isfile(matrix_file_full):
        logging.warning('Unable to find matrix file %s', matrix_file_full)
        return None
    cell_count = count_cells(matrix_file_full)
    logging.info('Cell count in %s is %s', accession_id, cell_count)
    return cell_count


def count_cells(matrix_file: str) -> Union[int, None]:
    """
    Count the number of cells in a matrix file

    :param matrix_file: A mtx file with tab/space/comma delimited data
    :return: A count of cells found in the given file
    """
    if not os.path.exists(matrix_file):
        logging.error(f'File not found "{matrix_file}"')
        return None
    first_line_length = len(open(matrix_file, newline='').readline().strip())
    if not first_line_length:
        logging.error(f'File has no first line "{matrix_file}"')
        return None
    barcode_indexes = set()
    with open(matrix_file, newline='') as csv_file:
        dialect = csv.Sniffer().sniff(csv_file.read(first_line_length))
        csv_file.seek(0)
        csv_reader = csv.reader(csv_file, dialect)
        next(csv_reader)  # skip header row
        for row in csv_reader:
            if len(row) == 3:
                barcode_indexes.add(row[1])  # 2nd column is barcode
    return len(barcode_indexes)


if __name__ == '__main__':
    main(sys.argv[1:])
