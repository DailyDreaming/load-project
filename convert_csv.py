import pandas as pd
import os


def csv_to_mtx(csv_filename: str, mtx_filename: str, cells_in_rows: bool) -> None:
    """
    :param csv_filename: path to input csv file.
    :param mtx_filename: path to output mtx directory.
    :param cells_in_rows: if true, cells are rows and genes are columns in the csv, else vice-versa.
    """
    csv = pd.read_csv(csv_filename, index_col=0)

    if not cells_in_rows:
        csv = csv.T

    cells = csv.index
    genes = csv.columns

    os.makedirs(mtx_filename, exist_ok=True)

    with open(os.path.join(mtx_filename, 'barcodes.tsv'), 'w') as f:
        f.writelines(cell + '\n' for cell in cells)

    with open(os.path.join(mtx_filename, 'genes.tsv'), 'w') as f:
        # CSV file sonly include gene names but scanpy requires a gene id/featurekey column so we fake it
        f.writelines('\t'.join([f'FAKE_FEATURE_KEY{i}', gene]) + '\n' for i, gene in enumerate(genes))

    entries = [
        (
            j + 1,
            i + 1,
            csv.iloc[i, j]
        )
        for i, _ in enumerate(cells)
        for j, _ in enumerate(genes)
        if csv.iloc[i, j] != 0
    ]

    def fmt_line(x):
        return ' '.join(map(str, x)) + '\n'

    with open(os.path.join(mtx_filename, 'matrix.mtx'), 'w') as f:
        f.write('%%MatrixMarket matrix coordinate real general\n')
        f.write(fmt_line(map(len, [genes, cells, entries])))
        f.writelines(map(fmt_line, entries))


def compile_mtxs(input_dir: str, output_filename: str) -> None:
    """
    :param input_dir: where to look for mtx files
    :param output_filename: path to output mtx directory
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

    merged = filetypes['genes.tsv'].merge(filetypes['cells.tsv']).merge(filetypes['matrix.mtx'])

    if merged.size != max(map(len, filetypes.values())):
        raise RuntimeError(f'There appear to be some missing/extra/unassociated mtx files')

    # merged now contains a single column of all the distinct matrix names (filename before extension)

    all_genes = set()

    def add_genes(prefix):
        with open(prefix + 'genes.tsv', 'r') as f:
            for line in f:
                all_genes.add(line.rstrip())

    merged.apply(add_genes,  axis=1)


if __name__ == '__main__':
    print('testing')

    import scanpy as sc
    import numpy as np

    data = np.random.random((100, 200))
    data[data < 0.5] = 0

    df = pd.DataFrame(data)
    df.columns = [str(chr(x)) + str(i) for x in range(ord('A'), ord('Z')) for i in range(10)][:200]
    df.index = [str(chr(x)) + str(i) for x in range(ord('a'), ord('z')) for i in range(5)][:100]

    df.to_csv('/tmp/test.csv')
    csv_to_mtx('/tmp/test.csv', '/tmp/test.mtx', True)

    adata1 = sc.read_csv('/tmp/test.csv')
    adata2 = sc.read_10x_mtx('/tmp/test.mtx')

    assert adata1.to_df().equals(adata2.to_df())
