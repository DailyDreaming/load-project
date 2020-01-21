from _pathlib import Path

from util import get_target_project_dirs


def populate_static_project(project_dir: Path, file_pattern: str):
    hca_dir = project_dir / 'hca'
    if hca_dir.is_dir():
        bundle_dir = project_dir / 'bundle'
        if bundle_dir.is_dir():
            for file in bundle_dir.glob(file_pattern):
                file.unlink()
        else:
            bundle_dir.mkdir()
        for file in hca_dir.glob(file_pattern):
            (bundle_dir / file.name).link_to(file)


def populate_all_static_projects(file_pattern):
    for project_dir in get_target_project_dirs():
        populate_static_project(project_dir, file_pattern)
