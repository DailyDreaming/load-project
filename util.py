import gzip
import os
from pathlib import Path
from typing import (
    Sequence,
    List,
    Optional,
)
import uuid


def open_maybe_gz(path, mode: str, **kwargs):
    """
    Like the open() built-in but transparently handles .gz files.

    Can be used as a context manager.
    """
    # Since `open` and `gzip.open` disagree on the default mode and whether just
    # 'r' implies 'b' or 't', the caller must be unambiguously specify the mode.
    # Plus, we don't support any of the write modes.
    if mode in ('rb', 'rt'):
        with open(path, 'rb') as f:
            header = f.read(2)
        if header == b'\x1f\x8b':  # GZIP magic cookie
            return gzip.open(path, mode, **kwargs)
        else:
            return open(path, mode, **kwargs)
    else:
        raise ValueError("Unsupported mode (must be 'rb' or 'rt'):", mode)


def generate_project_uuid(geo_accessions: Sequence[str]) -> str:
    """
    Deterministically generate a project UUID based on one or more GEO accession ids.
    """
    if isinstance(geo_accessions, str):
        geo_accessions = [geo_accessions]
    namespace_uuid = uuid.UUID('296e1002-1c99-4877-bb7f-bb6a3b752638')
    return str(uuid.uuid5(namespace_uuid, ''.join(sorted(geo_accessions))))


def generate_file_uuid(bundle_uuid: str, file_name: str) -> str:
    """
    Deterministically generate a file UUID based on the parent bundle uuid and its file name.
    """
    namespace_uuid = uuid.UUID('4c52e3d0-ffe5-4b4d-a4d0-cb6a6f372b31')
    return str(uuid.uuid5(namespace_uuid, bundle_uuid + file_name))


def get_skunk_accessions() -> Optional[List[str]]:
    try:
        accessions = os.environ['SKUNK_ACCESSIONS']
    except KeyError:
        return None
    else:
        return accessions.split(',')


def get_target_spreadsheets() -> List[Path]:
    accessions = get_skunk_accessions()
    spreadsheet_paths = []
    ext = '.0.xlsx'
    for sub_dir in ('existing', 'new'):
        src_dir = Path('spreadsheets') / sub_dir
        if accessions is None:
            subdir_paths = [p for p in src_dir.iterdir() if p.is_file() and p.name.endswith(ext)]
        else:
            subdir_paths = [src_dir / (accession + ext) for accession in accessions]
        spreadsheet_paths.extend(subdir_paths)
    return spreadsheet_paths


def get_target_project_dirs(uuids: bool = False, root_dir: Path = None) -> List[Path]:
    if root_dir is None:
        root_dir = Path('projects/')

    accessions = get_skunk_accessions()
    if accessions is None:
        return [p for p in root_dir.iterdir() if p.is_dir() and (p.is_symlink() ^ uuids)]
    else:
        targets = [generate_project_uuid([acc]) for acc in accessions] if uuids else accessions
        return [root_dir / target for target in targets]
