import argparse
import csv
import os
import sys

from typing import Sequence

import logging
logging.basicConfig(level=logging.INFO)


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--list', '-l',
                       metavar='MATRIX_FILE',
                       help='list cell count value for given matrix file')
    group.add_argument('--write', '-w',
                       metavar='MATRIX_FILE',
                       help='write a cell count file for given matrix file')
    group.add_argument('--list-all', '-L',
                       action='store_true',
                       help='list cell counts for all project matrix files')
    group.add_argument('--write-all', '-W',
                       action='store_true',
                       help='write cell count files for all project matrix files')
    args = parser.parse_args(argv)

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

    :param save_to_file: If true write each cell count to a file in the project's folder
    """
    for project_path in get_project_paths():
        logging.info('Checking: %s ...', project_path)
        count_one_project_cells(project_path, save_to_file)


def count_one_project_cells(path: str, save_to_file=True):
    """
    Count cells in one project

    :param path: A path to a project folder or matrix file
    :param save_to_file: If true write the cell count to a file in the project's folder
    """
    cell_count = get_project_cell_count(path)
    if save_to_file:
        if os.path.isdir(path):
            project_path = path
        else:
            project_path = '/'.join(path.split('/')[:-1])
        write_project_cell_count(project_path, cell_count)


def get_project_paths() -> Sequence[str]:
    """
    Return a list of the project folder paths
    """
    project_paths = []
    for item in os.listdir('projects'):
        project_path = f'projects/{item}'
        logging.debug('Checking: %s', project_path)
        if os.path.isdir(project_path) and not os.path.islink(project_path):
            project_paths.append(project_path)
    return project_paths


def write_project_cell_count(project_path: str, cell_count: int):
    """
    Save the cell count to a file in the project path

    :param project_path: A path to the project folder
    :param cell_count: A count of cells
    """
    if cell_count is not None:
        cell_count_file = f'{project_path}/cell_counts'
        logging.info('Writing: %s', cell_count_file)
        with open(cell_count_file, 'w') as f:
            f.write(f'{cell_count}\n')


def get_project_cell_count(path: str) -> int:
    """
    Get the cell count from a project's matrix file

    :param path: A path to a project folder or matrix file
    :return: A count of cells
    """
    matrix_file_full = None
    if os.path.isfile(path):
        matrix_file_full = path
    elif os.path.isdir(path):
        matrix_file_full = f'{path}/matrix.mtx'
    if not matrix_file_full or not os.path.exists(matrix_file_full):
        logging.debug('Unable to find matrix file at path: %s', path)
        return None
    cell_count = count_cells(matrix_file_full)
    logging.info('Cell count in %s: %s', path, cell_count)
    return cell_count


def count_cells(matrix_file: str) -> int:
    """
    Count the number of cells in a matrix file

    :param matrix_file: A mtx file with tab/space/comma delimited data
    :return: A count of cells found in the given file
    """
    if not os.path.exists(matrix_file):
        logging.error(f'File not found "{matrix_file}"')
        return None
    barcode_indexes = set(())
    with open(matrix_file, newline='') as csv_file:
        dialect = csv.Sniffer().sniff(csv_file.read(1024))
        csv_file.seek(0)
        csv_reader = csv.reader(csv_file, dialect)
        next(csv_reader)  # skip header row
        for row in csv_reader:
            if len(row) == 3:
                barcode_indexes.add(row[1])  # 2nd column is barcode
    return len(barcode_indexes)


if __name__ == '__main__':
    main(sys.argv[1:])
