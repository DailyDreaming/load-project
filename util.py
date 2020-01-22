import gzip
import os
from typing import (
    List,
    Optional,
    Sequence,
)
import uuid

from _pathlib import Path


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


WORKING_SET_ENV_VAR = 'SKUNK_ACCESSIONS'


def is_working_set_defined() -> bool:
    return WORKING_SET_ENV_VAR in os.environ


def get_skunk_accessions() -> Optional[List[str]]:
    try:
        accessions = os.environ[WORKING_SET_ENV_VAR]
    except KeyError:
        return None
    else:
        return [acc.strip() for acc in accessions.split(',') if acc.strip()]


def get_target_spreadsheets() -> List[Path]:
    accessions = get_skunk_accessions()
    spreadsheet_paths = []
    ext = '.0.xlsx'
    for sub_dir in ('existing', 'new'):
        src_dir = Path('spreadsheets') / sub_dir
        subdir_paths = [
            path
            for path
            in src_dir.iterdir()
            if (path.is_file()
                and path.name.endswith(ext)
                and (accessions is None
                     or path.name.replace(ext, '') in accessions))
        ]
        spreadsheet_paths.extend(subdir_paths)
    return spreadsheet_paths


def get_target_project_dirs(follow_links: bool = False) -> List[Path]:
    """
    Return all or a subset of the project directories, if that subset is
    configured.

    :param follow_links: If True, follow the symbolic accession links and return
                         Path instances referring to the physical, UUID-based
                         project directories. Otherwise Path instances referring
                         to the symbolic accession links will be returned.
    """
    projects_dir = Path('projects')
    accessions = get_skunk_accessions()
    symlinks = [
        path for path in projects_dir.iterdir()
        if path.is_dir() and path.is_symlink() and (accessions is None or path.name in accessions)
    ]
    if follow_links:
        project_dirs = []
        for symlink in symlinks:
            project_dir = symlink.follow()
            assert project_dir.is_dir() and not project_dir.is_symlink()
            accession = symlink.name
            project_uuid = generate_project_uuid([accession])
            assert project_dir.name == project_uuid
            project_dirs.append(project_dir)
        return project_dirs
    else:
        return symlinks
