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
    Iterable,
)
from contextlib import contextmanager
from pathlib import Path

import os
from more_itertools import one

log = logging.getLogger(__file__)


def synthesize_matrices(projects: Path):
    # GEO accessions for which the script will fail without special instructions
    # FIXME: These are not actual special cases, but are added for testing
    special_cases = {
        '06a318d9-54d8-5e41-aab5-f2d682fba690',
        '06f8848d-9c54-5829-92d3-d334809ad1e2',
        '096b7311-2bf7-5e61-9afb-d65c24a71243',
    }

    failed_projects = {}
    succeeded_projects = set()
    for project_dir in projects.iterdir():
        # Assuming that dir.name is the project UUID
        project_dir = Path(project_dir)
        project_uuid = extract_uuid(project_dir)
        if final_matrix_file(project_dir).exists():
            log.info('Final matrix already exists for project %s; moving on.', project_dir)
        elif project_uuid in special_cases:
            # FIXME
            log.warning('Skipping special case %s', project_dir)
            continue
        else:
            with tempfile.TemporaryDirectory() as staging_dir:
                try:
                    convert_and_move_matrices(project_dir, Path(staging_dir))
                except Exception as e:
                    failed_projects[project_uuid] = e
                    log.exception('Failed to process project', exc_info=True)
                else:
                    export_matrices(project_dir, staging_dir)
                    succeeded_projects.add(project_dir)

    print('Failed projects', file=sys.stderr)
    for p in failed_projects:
        print(p, file=sys.stderr)

    print('Succeeded projects')
    for project_dir in succeeded_projects:
        print(project_dir)


def convert_and_move_matrices(project_dir: Path, staging_dir: Path):
    # For the first iteration, just link everything that's in the geo dir into
    # the staging dir
    geo = geo_dir(project_dir)
    for dir_path, _, file_names in os.walk(geo):
        dir_path = Path(dir_path)
        for file in file_names:
            full_path = dir_path / file
            relative_path = full_path.relative_to(geo)
            os.makedirs(staging_dir / relative_path.parent, exist_ok=True)
            os.link(full_path, staging_dir / relative_path)


def find_mtx_files(filepaths: Iterable[Path]) -> Dict[str, Tuple[str, str, str]]:
    result = {}

    filenames = list(map(str, filepaths))

    anchors = [
        (fn, strip_suffix(strip_suffix(strip_suffix(fn, '.gz'), '.mtx'), 'matrix'))
        for fn, fp in zip(filenames, filepaths)
        if is_mtx(fp)
    ]

    for anchor_file, prefix in anchors:
        links = [
            file
            for file in filenames
            if file.startswith(prefix) and file != anchor_file
        ]
        if len(links) >= 2:
            def find(names):
                return one(link
                           for link in links
                           if any(strip_prefix(link, prefix).startswith(name) for name in names)
                           and strip_suffix(link, '.gz')[-4:] in ['.csv', '.tsv'])

            try:
                barcodes_file = find(['barcodes', 'cells'])
                genes_file = find(['genes', 'features'])
            except ValueError:
                log.warning('Could not identify row and column files for mtx')
            else:
                result[prefix] = (genes_file, barcodes_file, anchor_file)
        else:
            log.warning('There appear to be some missing/extra/unassociated mtx files')

    return result


def export_matrices(project_dir: Path, staging_dir: Path):
    """
    Zip up the final, converted matrices (in staging_dir) and copy the zip to
    its expected location.
    """
    final_matrix = final_matrix_file(project_dir)
    os.makedirs(final_matrix.parent, exist_ok=True)
    # make_archive adds it's own .zip at the end
    shutil.make_archive(strip_suffix(final_matrix_file(project_dir), '.zip'), 'zip', staging_dir)


def strip_suffix(s, suffix):
    s = str(s)
    if s.endswith(suffix):
        return s[:-len(suffix)]
    else:
        return s


def strip_prefix(s, prefix):
    s = str(s)
    if s.startswith(prefix):
        return s[len(prefix):]
    else:
        return s


def files_recursively(path) -> Sequence[Path]:
    for dir_path, _, files in os.walk(path):
        for f in files:
            yield Path(dir_path, f)


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


def geo_dir(project_dir: Path) -> Path:
    return project_dir / 'geo'


def bundle_dir(project_dir: Path) -> Path:
    return project_dir / 'bundle'


def final_matrix_file(project_dir: Path) -> Path:
    return bundle_dir(project_dir) / 'matrix.mtx.zip'


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


if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s',
                        level=logging.DEBUG)

    # To download some projects for testing run:
    # scp -r ubuntu@skunk.dev.explore.data.humancellatlas.org:/home/ubuntu/load-project/projects/0* ./test/projects
    synthesize_matrices(Path('test/projects'))
