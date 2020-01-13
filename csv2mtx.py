import argparse
import csv
import gzip
import io
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
)

import pandas as pd
from util import open_maybe_gz

logging.basicConfig(level=logging.INFO)

RowFilter = Callable[[List[str]], Optional[bool]]


class RowConverter(Iterable):

    def __init__(self,
                 row_provider: Iterable[List[str]],
                 rows_are_genes: bool,
                 row_filter: Optional[RowFilter] = None):
        """
        :param rows_are_genes: True if csv rows are genes, False if csv cols are barcodes
        :param row_filter: callable to process or skip nonconforming rows.
        """
        self.row_provider = row_provider
        self.rows_are_genes = rows_are_genes
        if row_filter is None:
            def row_filter(_): return False
        self.row_filter = row_filter
        self.x_axis_values = None
        self.y_axis_values = []
        self.num_values = 0

    def __iter__(self):

        for row in self.row_provider:
            filter_status = self.row_filter(row)
            if filter_status is None or filter_status is False:
                if self.x_axis_values is None:  # Get header values once
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
            elif filter_status is True:
                pass
            else:
                assert False, f"Invalid row_filter return type {type(filter_status)}"

    @property
    def genes(self):
        return self.y_axis_values if self.rows_are_genes else self.x_axis_values

    @property
    def barcodes(self):
        return self.x_axis_values if self.rows_are_genes else self.y_axis_values


class CSV2MTXConverter(RowConverter):
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
        """

        # Delay provider initialization until file is opened
        # noinspection PyTypeChecker
        super().__init__(None, rows_are_genes, row_filter)
        self.input_file = input_file
        self.delimiter = delimiter

    def __iter__(self):
        with open_maybe_gz(self.input_file, 'rt', newline='') as csv_file:
            self.row_provider = csv.reader(csv_file, delimiter=self.delimiter)
            return super().__iter__()


def convert_to_mtx(converter: RowConverter,
                   output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True) # FIXME: move to convert_matrices.py

    mtx_body_file = output_dir / 'matrix.mtx.body.gz'
    mtx_file = output_dir / 'matrix.mtx.gz'

    # Fully consume the iterator by writing the body of the mtx file to a temp file
    write_gzip_file(mtx_body_file, converter)

    # Write the completed mtx file using correct header information and the body we wrote to the temp file
    rows_cols_count_line = f'{len(converter.genes)} {len(converter.barcodes)} {converter.num_values}'
    write_mtx_file(rows_cols_count_line, mtx_body_file, mtx_file)
    mtx_body_file.unlink()

    # Write the two remaining files using the properties from the fully consumed iterator
    write_gzip_file(output_dir / 'barcodes.tsv.gz', ['barcodes'] + converter.barcodes)
    write_gzip_file(output_dir / 'genes.tsv.gz', ['genes'] + converter.genes)

    print('Done.')


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
    convert_to_mtx(csv_converter, output_dir)


def convert_cell_files_to_mtx(input_dir: Path,
                              output_dir: Path,
                              delimiter: str = ','):
    """
    Convert a series of files to mtx where each file is expression for a single cell.
    :param delimiter: what separates genes from their expression in the cell files.
    :param input_dir: directory containing the cell files.
    :param output_dir: The location to save the output files
    :return:
    """
    def iter_files():
        first = True
        for file_path in input_dir.iterdir():
            cell = pd.read_csv(file_path, sep=delimiter, compression='infer')
            if first:
                first = False
                # provide header (gene names) with empty first column
                genes = ['']
                genes.extend(cell[0])
                yield genes
            # provide expression values with barcodes
            data = [file_path.name]
            data.extend(cell[1])
            yield data

    cell_file_converter = RowConverter(iter_files(),
                                       rows_are_genes=False)

    convert_to_mtx(cell_file_converter, output_dir)


def write_gzip_file(output_file: Path, lines: Iterable):
    """
    Create/overwrite a gzipped text file

    :param output_file: File to create
    :param lines: List/Iterator of strings to write to file (a '\n' is added to each line)
    """
    temp_output_file = output_file.with_suffix(output_file.suffix + '.tmp')
    print(f'Writing {temp_output_file} ...')
    try:
        # Using gzip.open(temp) directly creates an archive that causes
        # `gunzip -N` to extract the file under the name of the temporary file
        # even if the archive name is different. Therefore we must set the
        # internal file name manually and pass in an already open file object
        # for writing.
        with open(str(temp_output_file), 'wb') as f:
            with gzip.GzipFile(filename=output_file, fileobj=f) as z:
                with io.TextIOWrapper(z) as w:
                    for line in lines:
                        w.write(line + '\n')
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
                # Using 1MiB buffer should be faster than the default of 16KiB
                copyfileobj(temp_data, f, length=2 ** 20)
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
