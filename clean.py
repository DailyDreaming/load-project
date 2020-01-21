import logging
import os
from _pathlib import Path
import shutil
from util import get_target_project_dirs

log = logging.getLogger(__name__)


def clean_project(project_dir: Path):
    log.warning(f'Removing matrix files for project {project_dir.name}')

    shutil.rmtree(project_dir / 'matrices', ignore_errors=True)
    try:
        os.remove(str(project_dir / 'bundle' / 'matrix.mtx.zip'))
    except FileNotFoundError:
        pass


def main():
    projects = get_target_project_dirs()
    for project in projects:
        clean_project(project)


if __name__ == "__main__":
    main()
