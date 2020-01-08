import argparse
import csv
import gzip
import logging
import os
from pathlib import Path
from shutil import copyfileobj
import sys
from typing import (
    Callable,
    Iterable,
    List,
    Optional,
    Union,
)

from util import open_maybe_gz

logging.basicConfig(level=logging.INFO)

RowFilter = Callable[[List[str]], None]


class CSV2MTXConverter(Iterable):
    """
    Convert a csv file to a matrix.mtx, barcodes.tsv, and genes.tsv set of files
    """

    def __init__(self,
                 input_file: Path,
                 delimiter: str = ',',
                 rows_are_genes: bool = True,
                 row_filter: Optional[RowFilter] = None):
        """
        :param input_file: The input csv file
        :param delimiter: Delimiter character in csv
        :param rows_are_genes: True if csv rows are genes, False if csv cols are barcodes
        """
        self.input_file = input_file
        self.delimiter = delimiter
        self.rows_are_genes = rows_are_genes
        if row_filter is None:
            def row_filter(_):
                pass
        self.row_filter = row_filter
        self.x_axis_values = []
        self.y_axis_values = []
        self.num_values = 0

    def __iter__(self):
        with open_maybe_gz(self.input_file, 'rt', newline='') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=self.delimiter)
            for row in csv_reader:
                filter_status = self.row_filter(row)
                if not filter_status:  # A True filter_status tells us to skip the entire row
                    if not self.x_axis_values:  # Get header values once
                        self.x_axis_values = row[1:]
                    else:
                        self.y_axis_values.append(row[0])
                        for col, value in enumerate(row[1:]):
                            if float(value):
                                self.num_values += 1
                                gene_index = len(self.y_axis_values) if self.rows_are_genes else col + 1
                                barcode_index = col + 1 if self.rows_are_genes else len(self.y_axis_values)
                                return_string = f'{gene_index} {barcode_index} {value}'
                                yield return_string

    @property
    def genes(self):
        return self.y_axis_values if self.rows_are_genes else self.x_axis_values

    @property
    def barcodes(self):
        return self.x_axis_values if self.rows_are_genes else self.y_axis_values


def convert_csv_to_mtx(input_file: Path,
                       output_dir: Path,
                       delimiter: str = ',',
                       rows_are_genes: bool = True,
                       row_filter: Optional[RowFilter] = None):
    """
    Convert a csv file to a matrix.mtx, barcodes.tsv, and genes.tsv set of files

    :param input_file: The input csv file
    :param output_dir: The location to save the output files
    :param delimiter: Delimiter character in csv
    :param rows_are_genes: True if csv rows are genes, False if csv cols are barcodes
    :param row_filter: callable to pre-process lines in the CSV
    """
    csv_converter = CSV2MTXConverter(input_file, delimiter, rows_are_genes, row_filter)
    output_dir.mkdir(parents=True, exist_ok=True)

    mtx_body_file = output_dir / 'matrix.mtx.body.gz'
    mtx_file = output_dir / 'matrix.mtx.gz'

    # Fully consume the iterator by writing the body of the mtx file to a temp file
    write_gzip_file(mtx_body_file, csv_converter)

    # Write the completed mtx file using correct header information and the body we wrote to the temp file
    rows_cols_count_line = f'{len(csv_converter.genes)} {len(csv_converter.barcodes)} {csv_converter.num_values}'
    write_mtx_file(rows_cols_count_line, mtx_body_file, mtx_file)
    mtx_body_file.unlink()

    # Write the two remaining files using the properties from the fully consumed iterator
    write_gzip_file(output_dir / 'barcodes.tsv.gz', ['barcodes'] + csv_converter.barcodes)
    write_gzip_file(output_dir / 'genes.tsv.gz', ['genes'] + csv_converter.genes)

    print('Done.')


def write_gzip_file(output_file: Path, lines: Union[Iterable, list]):
    """
    Create/overwrite a gzipped text file

    :param output_file: File to create
    :param lines: List/Iterator of strings to write to file (a '\n' is added to each line)
    """
    temp_output_file = output_file.with_suffix(output_file.suffix + '.tmp')
    print(f'Writing {temp_output_file} ...')
    try:
        with gzip.open(temp_output_file, 'wb') as f:
            for line in lines:
                f.write((str(line) + '\n').encode())
    except:
        try:
            temp_output_file.unlink()
        except FileNotFoundError:
            pass
        raise
    else:
        print(f'Renaming {temp_output_file} to {output_file} ...')
        temp_output_file.rename(output_file)


def write_mtx_file(rows_cols_count_line: str, mtx_body_file: Path, output_file: Path):
    """
    Write the final mtx file with comment header line, the rows_cols_count line, and
    the mtx body from previously written temp file

    :param rows_cols_count_line: String containing "{num_genes} {num_cells} {total_values}"
    :param mtx_body_file: Path of the temp file containing data to be written to the body mtx file
    :param output_file: Path of the mtx file to be written
    """
    temp_output_file = output_file.with_suffix(output_file.suffix + '.tmp')
    try:
        with gzip.open(temp_output_file, 'wb') as f:
            header_line = '%%MatrixMarket matrix coordinate integer general\n'
            f.write(header_line.encode())
            f.write((rows_cols_count_line + '\n').encode())
            with open_maybe_gz(mtx_body_file, 'rb') as temp_data:
                copyfileobj(temp_data, f)
    except:
        try:
            temp_output_file.unlink()
        except FileNotFoundError:
            pass
        raise
    else:
        print(f'Renaming {temp_output_file} to {output_file} ...')
        temp_output_file.rename(output_file)


def main(argv):
    """
    Support for command line execution of convert_csv_to_mtx()
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('csv_file', help='Input csv file')
    parser.add_argument('output_dir', help='Path to write output files')
    parser.add_argument('delimiter', help='Delimiter character or keyword "comma", "space", "tab"')
    parser.add_argument('rows_are_genes', help='"y" if rows are genes or "n" if columns are genes')
    args = parser.parse_args(argv)

    if not os.path.isfile(args.csv_file):
        print('Error: File not found:', args.csv_file)
        parser.print_help()
        exit()

    if args.delimiter == 'comma':
        args.delimiter = ','
    elif args.delimiter == 'space':
        args.delimiter = ' '
    elif args.delimiter == 'tab':
        args.delimiter = '\t'

    if len(args.delimiter) < 1:
        print('Error: delimiter must be 1 char in length')

    if args.rows_are_genes not in ('y', 'n'):
        print('Error: rows_are_genes must be "y" or "n"')

    args.rows_are_genes = args.rows_are_genes == 'y'

    convert_csv_to_mtx(Path(args.csv_file),
                       Path(args.output_dir),
                       delimiter=args.delimiter,
                       rows_are_genes=args.rows_are_genes)


if __name__ == '__main__':
    main(sys.argv[1:])
