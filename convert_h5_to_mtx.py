import sys
import os
import argparse
import h5py
import time


def simple_validation(count, file):
    with open(file, 'r') as f:
        lines = [l for l in f]
    assert len(lines) == count, f'{file} {len(lines)} != expected {count}'


def run(input_hd5, output_dir):
    output_mtx = os.path.join(output_dir, 'matrix.mtx')
    output_barcodes = os.path.join(output_dir, 'barcodes.tsv')
    output_genes = os.path.join(output_dir, 'genes.tsv')

    overall_start = time.time()

    with h5py.File(input_hd5, 'r') as h5_file, \
            open(output_mtx, 'w') as mtx_file, \
            open(output_barcodes, 'w') as barcodes_file, \
            open(output_genes, 'w') as genes_file:
        assert len(list(h5_file.keys())) == 1
        group = list(h5_file.keys())[0]

        # ['barcodes', 'data', 'gene_names', 'genes', 'indices', 'indptr', 'shape']
        # https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices
        data = list(h5_file[group])
        assert sorted(data) == sorted(['barcodes', 'data', 'gene_names', 'genes', 'indices', 'indptr', 'shape'])

        number_of_genes = h5_file[group]['gene_names'].shape[0]
        number_of_cells = h5_file[group]['barcodes'].shape[0]
        number_of_total_data_points = h5_file[group]['data'].shape[0]

        # write the genes file
        start = time.time()
        for i, gene in enumerate(h5_file[group]['genes']):
            gene = gene.decode("utf-8")
            gene_name = h5_file[group]['gene_names'][i].decode("utf-8")
            genes_file.write(f'{gene}\t{gene_name}\n')
        diff = time.time() - start
        print(f'Genes File Complete.  Time Taken: {diff} seconds.')

        # write the barcodes file
        start = time.time()
        for barcode in h5_file[group]['barcodes']:
            barcodes_file.write(f'{barcode.decode("utf-8")}\n')
        diff = time.time() - start
        print(f'Barcodes File Complete.  Time Taken: {diff} seconds.')

        # write the main mtx file
        mtx_header = f'%%MatrixMarket matrix coordinate integer general\n' \
                     f'%\n' \
                     f'{number_of_genes} {number_of_cells} {number_of_total_data_points}\n'
        mtx_file.write(mtx_header)

        sample_range = iter(h5_file[group]['indptr'][:])
        current_range = next(sample_range)
        cell_barcode_id = 0
        start = time.time()
        for i, gene_expression_level in enumerate(h5_file[group]['data']):
            if i == current_range and i != number_of_total_data_points:
                current_range = next(sample_range)
                cell_barcode_id += 1

            demarcation = int(number_of_total_data_points / 100) + 1
            if (i % demarcation == 0 or i == number_of_total_data_points) and i != 0:
                diff = time.time() - start
                time_taken = int(diff / 6) / 10
                time_required = int((diff / (i / number_of_total_data_points)) / 6) / 10
                print(f'{int(i * 100 / number_of_total_data_points)}% Completed.  Time taken: {time_taken} minutes / {time_required} estimated minutes total.')
            gene_id = h5_file[group]['indices'][i]
            gene_coord = f'{gene_id} {cell_barcode_id} {gene_expression_level}\n'
            mtx_file.write(gene_coord)

    diff = time.time() - overall_start
    time_taken = int(diff / 6) / 10
    print(f'Complete.  Total Time Taken: {time_taken} minutes.')

    # sanity checks
    simple_validation(count=number_of_cells, file=output_barcodes)
    simple_validation(count=number_of_genes, file=output_genes)
    header_lines = len(mtx_header.strip().split('\n'))
    simple_validation(count=number_of_total_data_points + header_lines, file=output_mtx)


def main(argv):
    parser = argparse.ArgumentParser(description='md5 format to mtx format.')
    parser.add_argument("--hd5", type=str,
                        default='projects/458cbaeb-c4b7-5537-b1f3-a5d537478112/geo/GSE118127_RAW/GSM3319032_sample_1-1_filtered_gene_bc_matrices_h5.h5',
                        help="Path to an hd5 file.  Example: 'data/test.h5'")
    parser.add_argument("--output_dir", type=str,
                        default='/home/quokka/Desktop/mtx/output/',
                        help="Path to an output directory.")

    args = parser.parse_args(argv)
    print(f'Beginning conversion of: {args.hd5}')
    run(input_hd5=args.hd5, output_dir=args.output_dir)


if __name__ == "__main__":
    main(sys.argv[1:])
