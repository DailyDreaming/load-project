import collections
import gzip
import logging
import mimetypes
import shutil
import sys
import tempfile
from typing import (
    Tuple,
    Dict,
    Sequence,
)
import zipfile
from contextlib import contextmanager
from pathlib import Path

import pandas as pd
import numpy as np
import os
from more_itertools import one

import conversions

log = logging.getLogger(__file__)


def try_delim_file(file: str, target: str, hints: dict = None) -> None:
    hints = {} if hints is None else hints
    try:
        # Pandas will infer sep if None
        data = pd.read_csv(file, index_col=0, sep=hints.get('separator'))
    except Exception:
        log.warning(f'{file}: not a parse-able data file')
        raise

    if data.shape[1] == 1:
        log.warning(f'{file}: not a matrix file, or using wrong separator.')
        raise RuntimeError

    cells_in_rows = hints.get('cells_in_rows')

    # Try to guess which axis is genes based on size.
    # Humans and mice both have approx. 25000 genes so one axis is close to this size and the other is far from it we
    # can fairly certain which is genes
    if cells_in_rows is None:
        n_genes = 25000
        genes_in_rows = abs(data.shape[0] - n_genes) <= 5000
        cells_in_rows = abs(data.shape[1] - n_genes) <= 5000

        # other possible strategies for heuristic:
        # search for a gene name in index or columns

        if cells_in_rows == genes_in_rows:
            log.warning(f'{file}: could not infer which axis was genes (size = {data.shape})')
            raise RuntimeError
        else:
            log.info(f'{file}: inferred cells are in {"rows" if cells_in_rows else "columns"} from shape {data.shape}')

    if cells_in_rows is False:
        data = data.T

    try:
        csv_to_mtx(data, target)
    except Exception:
        log.warning(f'{file}: conversion step failed')
        raise


def csv_to_mtx(data: pd.DataFrame, mtx_filename: str) -> None:
    cells = data.index
    genes = data.columns

    # To properly cross-reference this we'd need to know species and maybe other info
    if all(gene.startswith('ENS') for gene in genes[10:]):
        gene_ids = genes
        gene_names = [f'FAKE_GENE_NAME_{gene}' for gene in genes]
    else:
        gene_ids = [f'FAKE_GENE_ID_{gene}' for gene in genes]
        gene_names = genes

    os.makedirs(mtx_filename, exist_ok=True)

    with open(os.path.join(mtx_filename, 'barcodes.tsv'), 'w') as f:
        f.writelines(cell + '\n' for cell in cells)

    with open(os.path.join(mtx_filename, 'genes.tsv'), 'w') as f:
        # CSV files only include gene names but scanpy requires a gene id/featurekey column so we fake it
        # Use hash for consistency across files
        f.writelines('\t'.join([gene_id, gene_name]) + '\n' for gene_id, gene_name in zip(gene_ids, gene_names))

    entries = [
        (
            gene_idx + 1,
            cell_idx + 1,
            data.iloc[cell_idx, gene_idx]
        )
        for cell_idx, _ in enumerate(cells)
        for gene_idx, _ in enumerate(genes)
        if data.iloc[cell_idx, gene_idx] != 0
    ]

    def fmt_line(x):
        return ' '.join(map(str, x)) + '\n'

    with open(os.path.join(mtx_filename, 'matrix.mtx'), 'w') as f:
        f.write('%%MatrixMarket matrix coordinate real general\n')
        f.write(fmt_line(map(len, [genes, cells, entries])))
        f.writelines(map(fmt_line, entries))


def find_mtx_files(input_dir: Path) -> Dict[str, Tuple[str, str, str]]:
    result = {}

    filepaths = input_dir.iterdir()
    filenames = list(map(str, filepaths))

    anchors = [
        (fn, strip_suffix(strip_suffix(strip_suffix(fn, '.gz'), '.mtx'), 'matrix'))
        for fn, fp in zip(filenames, filepaths)
        if is_mtx(fp)
    ]

    for anchor_file, prefix in anchors:
        links = [
            strip_prefix(file, prefix)
            for file in filenames
            if file.startswith(prefix) and file != anchor_file
        ]
        if len(links) >= 2:
            def find(names):
                return one(link
                           for link in links
                           if any(link.startswith(name) for name in names)
                           and strip_suffix(link, '.gz')[-4:] in ['.csv', '.tsv'])

            try:
                barcodes_file = find(['barcodes', 'cells'])
                genes_file = find(['genes', 'features'])
            except ValueError:
                log.warning('Couldn\'t identify row and column files for mtx')
            else:
                result[prefix] = (genes_file, barcodes_file, anchor_file)
        else:
            log.warning('There appear to be some missing/extra/unassociated mtx files')

    return result


def compile_mtxs(triplets: Sequence[Tuple[str, str, str]], output_filename: str, tsv_headers: bool) -> None:
    """
    :param triplets: sequence of genes, barcodes, matrix files.
    :param output_filename: path to output mtx directory
    :param tsv_headers: whether there are column headers in the genes and barcodes files
    """

    entries = []

    def add_matrix(triplet):

        mtx = pd.read_csv(triplet[2], skiprows=2, sep=' ', header=None)
        mtx.columns = ['gene_idx', 'cell_idx', 'value']

        cells = pd.read_csv(triplet[1],
                            sep='\t',  # some of these are named .csv but I'm pretty sure they all actually use tabs
                            header=0 if tsv_headers else None)
        cells = cells.iloc[:, :1]  # only barcodes column
        cells.columns = ['barcode']
        cells.reset_index(inplace=True)
        cells['index'] += 1

        genes = pd.read_csv(triplet[0],
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
    map(add_matrix, triplets)
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
    genes.to_csv(os.path.join(output_filename, 'genes.tsv.gz'),
                 compression='gzip',
                 header=False,
                 index=False,
                 sep='\t')

    barcodes = entries[['barcode']].drop_duplicates()
    barcodes.to_csv(os.path.join(output_filename, 'barcodes.tsv.gz'),
                    compression='gzip',
                    header=False,
                    index=False,
                    sep='\t')

    entries = entries[['gene_idx', 'cell_idx', 'value']]
    mtx = os.path.join(output_filename, 'matrix.mtx')

    # appending to gzip archive doesn't seem to work
    with open(mtx, '+') as f:
        f.write('%%MatrixMarket matrix coordinate real general\n')
        f.write(' '.join([str(df.shape[0]) for df in (genes, barcodes, entries)]) + '\n')
        entries.to_csv(f, mode='a', header=False, index=False, sep=' ')

        with gzip.open(mtx + '.gz', 'wb') as gzf:
            shutil.copyfileobj(f, gzf)


def extract_uuid(path: Path):
    """
    Gets the project UUID from the path

    >>> extract_uuid(Path('~/load-project.jesse/projects/51a21599-a014-5c5a-9760-d5bdeb80f741/geo'))
    51a21599-a014-5c5a-9760-d5bdeb80f741
    """
    uuid_index = path.parts.index('projects') + 1
    return path.parts[uuid_index]


def is_mtx(p: Path):
    return p.name.endswith('.mtx') or p.name.endswith('.mtx.gz')


def geo_dir(project_dir: Path):
    return project_dir / 'geo'


def matrix_dir(project_dir: Path):
    return project_dir / 'matrix'


def staging_dir(project_dir: Path):
    return project_dir / 'matrix_staging'


def bundle_dir(project_dir: Path):
    return project_dir / 'bundle'


@contextmanager
def read_maybe_gz(filename, **kwargs):
    m_type, encoding = mimetypes.guess_type(filename)
    if encoding == 'gzip':
        open_ = gzip.open(filename, 'rt', encoding='utf-8', **kwargs)
    else:
        open_ = open(filename, 'r', newline='', **kwargs)
    try:
        with open_ as f:
            yield f
    except UnicodeDecodeError:
        log.warning(f'Cannot open `{filename}` since it is not text nor gzip. Maybe tar?')


def strip_suffix(s, suffix):
    if s.endswith(suffix):
        return s[:-len(suffix)]
    else:
        return s


def strip_prefix(s, prefix):
    if s.startswith(prefix):
        return s[len(prefix):]
    else:
        return s


def matrix_name(filename: Path):
    """
    This is more or less arbitrary. Just used to avoid namespace collisions from multiple
    converted matrices.
    """
    name = filename.name
    name = strip_suffix(name, '.gz')
    for suffix in ('.txt', '.csv', '.tsv'):
        name = strip_suffix(name, suffix)
    return name + '.mtx'


def files_recursively(path: Path) -> Sequence[Path]:
    for dir_path, _, files in os.walk(path):
        for f in files:
            yield Path(dir_path, f)


def try_to_convert(files: Sequence[Path], tmpdir: str) -> None:
    delim_exts = ['csv', 'tsv', 'txt']
    delim_exts.extend([ext + '.gz' for ext in delim_exts])

    for file in files:
        if any(file.name.endswith(ext) for ext in delim_exts):
            log.info(f'Found potential csv/tsv {file}, attempting to convert')
            try:
                try_delim_file(str(file),
                               Path(tmpdir) / matrix_name(Path(file)))
            except Exception:
                log.warning(f'Failed to convert file {file}', exc_info=True)


def synthesize_matrix(project_dir: Path):
    with tempfile.TemporaryDirectory() as tmpdir:
        project_uuid = extract_uuid(project_dir)
        if project_uuid == '099c02da-23b2-5748-8618-92bc6770dc51':
            # one mtx, just need to link files with proper names
            conversions.one_mtx(project_dir, ('GSE106273_combined_genes.tsv.gz',
                                              'GSE106273_combined_barcodes.tsv.gz',
                                              'GSE106273_combined_matrix.tsv.gz'))
        elif project_uuid == '0a8f2289-5862-5bf0-8c27-0885453de788':
            # multiple mtx files with "correct" names, no special case required
            default_synthesis_technique(project_dir, tmpdir)
        elif project_uuid == '061ec9d5-9acf-54db-9eee-555136d5ce41':
            # multiple mtxs with un-parseable filenames
            compile_mtxs([
                ('GSM3271040_RNA_sciCAR_A549_gene.txt.gz',  'GSM3271040_RNA_sciCAR_A549_cell.txt.gz',  'GSM3271040_RNA_sciCAR_A549_gene_count.txt.gz'),
                ('GSM3271042_RNA_only_A549_gene.txt.gz',    'GSM3271042_RNA_only_A549_cell.txt.gz',    'GSM3271042_RNA_only_A549_gene_count.txt.gz'),
                ('GSM3271044_RNA_mouse_kidney_gene.txt.gz', 'GSM3271044_RNA_mouse_kidney_cell.txt.gz', 'GSM3271044_RNA_mouse_kidney_gene_count.txt.gz')
            ], str(project_dir), True)
        else:
            # default_synthesis_technique(project_dir, tmpdir)
            log.info('Do the default thing')
            raise RuntimeError


def default_synthesis_technique(project_dir: Path, tmpdir: str):
    """
    Look at files in project and decide how they should be combined
    / transformed in order to produce a single .mtx file.

    The algorithm is:
    If we don't have any mtx files:
        Try and make some in a temp dir.
        Then squish them together.
    Otherwise:
        just squish together what we have.

    This assumes (probably correctly) that we won't encounter projects
    containing data in BOTH csv/tsv/txt AND mtx files.

    If unable to synthesize the matrices raises a RuntimeError
    """
    log.info('Starting project %s', project_dir)
    in_dir = geo_dir(project_dir)
    files = list(files_recursively(in_dir))
    log.info('Project %s contains %s files', project_dir, len(files))

    if any(is_mtx(f) for f in files):
        log.info('Found some mtx files; trying to convert them')
        mtxs = find_mtx_files(in_dir)
        # Incredibly, headers are only present in the CSV files, not the TSV ones.
        # THIS MAY BREAK ON NEW FILES!!
        tsv_headers = [strip_suffix(mtx[0], '.gz').endswith('.csv') for mtx in mtxs]
        if any(tsv_headers) and not all(tsv_headers):
            raise RuntimeError('Mixed csv and tsv files in mtx directory')
        compile_mtxs(list(mtxs.values()), str(matrix_dir(project_dir)), tsv_headers=tsv_headers[0])
    else:
        try_to_convert(files, tmpdir)
        mtxs = find_mtx_files(in_dir)
        files = list(files_recursively(tmpdir))
        if any(is_mtx(f) for f in files):
            log.info('%s files were converted to mtx; compiling', len(files))
            # Here we assume there are no TSV headers because csv_to_mtx doesn't add them.
            compile_mtxs(list(mtxs.values()), str(matrix_dir(project_dir)), tsv_headers=False)
        else:
            raise RuntimeError("Unable to synthesize; nothing converted / nothing to convert")


def zip_matrix(project_dir: Path):
    final_matrix = final_matrix_file(project_dir)
    os.makedirs(final_matrix.parent, exist_ok=True)
    with zipfile.ZipFile(str(final_matrix), 'w') as zipf:
        for filename in ['matrix.mtx.gz', 'barcodes.tsv.gz', 'genes.tsv.gz']:
            zipf.write(project_dir / filename, arcname=filename)


def final_matrix_file(project_dir: Path):
    return bundle_dir(project_dir) / 'matrix.mtx.zip'


def synthesize_matrices(projects: Path):
    # GEO accessions for which the script will fail without special instructions
    special_cases = {}

    failed_projects = {}
    succeeded_projects = set()
    for project_dir in projects.iterdir():
        # Assuming that dir.name is the project UUID
        project_dir = Path(project_dir)
        project_uuid = extract_uuid(project_dir)
        if not final_matrix_file(project_dir).exists():
            try:
                if project_uuid in special_cases:
                    special_cases[project_uuid](project_dir)
                else:
                    synthesize_matrix(project_dir)
            except Exception as e:
                failed_projects[project_uuid] = e
                log.exception('Failed to process project', exc_info=True)
            else:
                succeeded_projects.add(project_dir)

    print('Failed projects', file=sys.stderr)
    for p in failed_projects:
        print(p, file=sys.stderr)

    for project_dir in succeeded_projects:
        zip_matrix(project_dir)


if __name__ == '__main__':
    log.setLevel('INFO')
    # Noah, to download some test files run
    # scp -r ubuntu@skunk.dev.explore.data.humancellatlas.org:/home/ubuntu/load-project/projects/0* ./test/projects
    synthesize_matrices(Path('test/projects'))


def test_conversion():
    print('testing')

    import scanpy as sc
    import numpy as np

    data = np.random.random((100, 200))
    data[data < 0.5] = 0  # add a bit of sparseness

    df = pd.DataFrame(data)
    df.columns = [str(chr(x)) + str(i) for x in range(ord('A'), ord('Z')) for i in range(10)][:200]
    df.index = [str(chr(x)) + str(i) for x in range(ord('a'), ord('z')) for i in range(5)][:100]

    df.to_csv('/tmp/test.csv')
    csv_to_mtx('/tmp/test.csv', '/tmp/test.mtx')

    adata1 = sc.read_csv('/tmp/test.csv')
    adata2 = sc.read_10x_mtx('/tmp/test.mtx')

    assert adata1.to_df().equals(adata2.to_df())

    # some manually created files in here
    compile_mtxs('mtx_test_files', 'mtx_test_files/out/', False)
    adata = sc.read_10x_mtx('/tmp/mtxs/out')

    print(adata.to_df())

    # this doesn't really verify integrity but it does verify consistency
