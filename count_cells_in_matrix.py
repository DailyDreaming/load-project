import argparse
import csv
import os
import sys

import logging
logging.basicConfig(level=logging.INFO)


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('files',
                        nargs='+',
                        help='matrix.mtx file')
    args = parser.parse_args(argv)

    for file in args.files:
        logging.info(file + ": " + str(count_cells(file)) + ' cells')


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
