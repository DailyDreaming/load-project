import collections

import pandas as pd
import os
from more_itertools import one


def csv_to_mtx(csv_filename: str, mtx_filename: str, cells_in_rows: bool, **pd_args) -> None:
    """
    :param csv_filename: path to input csv file.
    :param mtx_filename: path to output mtx directory.
    :param cells_in_rows: if true, cells are rows and genes are columns in the csv, else vice-versa.
    :param pd_args: keywords args to be passed to pandas.read_csv
    """
    csv = pd.read_csv(csv_filename, index_col=0, **pd_args)

    if not cells_in_rows:
        csv = csv.T

    cells = csv.index
    genes = csv.columns

    os.makedirs(mtx_filename, exist_ok=True)

    with open(os.path.join(mtx_filename, 'barcodes.tsv'), 'w') as f:
        f.writelines(cell + '\n' for cell in cells)

    with open(os.path.join(mtx_filename, 'genes.tsv'), 'w') as f:
        # CSV file sonly include gene names but scanpy requires a gene id/featurekey column so we fake it
        # Use hash for consistency across files
        f.writelines('\t'.join([f'FAKE_FEATURE_KEY{hash(gene)}', gene]) + '\n' for gene in genes)

    entries = [
        (
            gene_idx + 1,
            cell_idx + 1,
            csv.iloc[cell_idx, gene_idx]
        )
        for cell_idx, _ in enumerate(cells)
        for gene_idx, _ in enumerate(genes)
        if csv.iloc[cell_idx, gene_idx] != 0
    ]

    def fmt_line(x):
        return ' '.join(map(str, x)) + '\n'

    with open(os.path.join(mtx_filename, 'matrix.mtx'), 'w') as f:
        f.write('%%MatrixMarket matrix coordinate real general\n')
        f.write(fmt_line(map(len, [genes, cells, entries])))
        f.writelines(map(fmt_line, entries))


def compile_mtxs(input_dir: str, output_filename: str, tsv_headers: bool) -> None:
    """
    :param input_dir: where to look for mtx files
    :param output_filename: path to output mtx directory
    :param tsv_headers: whether there are column headers in genes.tsv and barcodes.tsv
    """
    files = os.listdir(input_dir)
    fileexts = ('genes.tsv', 'barcodes.tsv', 'matrix.mtx')

    filetypes = {
        fileext: pd.DataFrame([
            [
                file.split(fileext)[0]
            ]
            for file in files
            if file.endswith(fileext)
        ])
        for fileext in fileexts
    }

    file_prefixes = filetypes['genes.tsv'].merge(filetypes['barcodes.tsv']).merge(filetypes['matrix.mtx'])

    if file_prefixes.size != max(map(len, filetypes.values())):
        # if this causes problems, an alternative would be to disregard unassociated files
        # and only compile ones where the full triplet is present.
        raise RuntimeError(f'There appear to be some missing/extra/unassociated mtx files')

    # file_prefixes now contains a single column of all the distinct matrix names (filename before extension)

    entries = []

    def add_matrix(row):
        prefix = one(row)

        mtx = pd.read_csv(os.path.join(input_dir, prefix + 'matrix.mtx'), skiprows=2, sep=' ', header=None)
        mtx.columns = ['gene_idx', 'cell_idx', 'value']

        cells = pd.read_csv(os.path.join(input_dir, prefix + 'barcodes.tsv'),
                            sep='\t',
                            header=0 if tsv_headers else None)
        cells = cells.iloc[:, :1]  # only barcodes column
        cells.columns = ['barcode']
        cells.reset_index(inplace=True)
        cells['index'] += 1

        genes = pd.read_csv(os.path.join(input_dir, prefix + 'genes.tsv'),
                            sep='\t',
                            header=0 if tsv_headers else None)
        genes = genes.iloc[:, :2]  # only id and name columns
        genes.columns = ['gene_id', 'gene_name']
        genes.reset_index(inplace=True)
        genes['index'] += 1

        mtx = mtx.merge(genes, left_on='gene_idx', right_on='index', validate='m:1')
        mtx.drop(['gene_idx', 'index'], inplace=True, axis=1)
        mtx = mtx.merge(cells, left_on='cell_idx', right_on='index', validate='m:1')
        mtx.drop(['cell_idx', 'index'], inplace=True, axis=1)

        entries.append(mtx)

    # stack matrix entries into one frame
    file_prefixes.apply(add_matrix, axis=1)
    entries = pd.concat(entries, axis=0, ignore_index=True)

    # consolidate inconsistent gene ids since most of them will probably be garbage from the previous method
    # this step might no longer be necessary given the hash changes
    for gene_name, group in entries.groupby('gene_name'):
        if group['gene_id'].nunique() > 1:
            mask = (entries['gene_name'] == gene_name)
            # just take the first one
            entries.loc[mask, 'gene_id'] = entries.loc[np.flatnonzero(mask)[0], 'gene_id']

    # resolve duplicate barcodes
    for barcode, group in entries.groupby('barcode'):
        gene_counts = group['gene_name'].value_counts()
        for gene in gene_counts.index:
            if gene_counts[gene] > 1:
                mask = (entries['barcode'] == barcode) & (entries['gene_name'] == gene)
                # take arithmetic mean of gene expression values
                entries.loc[mask, 'value'] = entries.loc[mask, 'value'].mean()

    # now that all duplicate rows have the same value we just drop the excess
    entries.drop_duplicates(subset=['barcode', 'gene_name'], inplace=True)

    # assign line numbers to unique gene names and barcodes
    class id_dict(collections.defaultdict):
        def __init__(self):
            super().__init__(lambda: 1 + len(self))

    entries['gene_idx'] = entries['gene_name'].map(id_dict())
    entries['cell_idx'] = entries['barcode'].map(id_dict())

    os.makedirs(output_filename, exist_ok=True)

    genes = entries[['gene_id', 'gene_name']].drop_duplicates()
    genes.to_csv(os.path.join(output_filename, 'genes.tsv'),
                 header=False,
                 index=False,
                 sep='\t')

    barcodes = entries[['barcode']].drop_duplicates()
    barcodes.to_csv(os.path.join(output_filename, 'barcodes.tsv'),
                    header=False,
                    index=False,
                    sep='\t')

    entries = entries[['gene_idx', 'cell_idx', 'value']]
    mtx = os.path.join(output_filename, 'matrix.mtx')
    with open(mtx, 'w') as f:
        f.write('%%MatrixMarket matrix coordinate real general\n')
        f.write(' '.join([str(df.shape[0]) for df in (genes, barcodes, entries)]) + '\n')
    entries.to_csv(mtx, mode='a', header=False, index=False, sep=' ')


if __name__ == '__main__':
    print('testing')

    import scanpy as sc
    import numpy as np

    data = np.random.random((100, 200))
    data[data < 0.5] = 0  # add a bit of sparseness

    df = pd.DataFrame(data)
    df.columns = [str(chr(x)) + str(i) for x in range(ord('A'), ord('Z')) for i in range(10)][:200]
    df.index = [str(chr(x)) + str(i) for x in range(ord('a'), ord('z')) for i in range(5)][:100]

    df.to_csv('/tmp/test.csv')
    csv_to_mtx('/tmp/test.csv', '/tmp/test.mtx', True)

    adata1 = sc.read_csv('/tmp/test.csv')
    adata2 = sc.read_10x_mtx('/tmp/test.mtx')

    assert adata1.to_df().equals(adata2.to_df())

    # some manullay created files in here
    compile_mtxs('/tmp/mtxs', '/tmp/mtxs/out/', False)
    adata = sc.read_10x_mtx('/tmp/mtxs/out')

    print(adata.to_df())

    # this doesn't really verify integrity but it does verify consistency
