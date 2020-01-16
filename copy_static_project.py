import os
import re
import glob
from pathlib import Path

UUID_PATTERN = "[a-z0-9]{8}-[a-z0-9]{4}-[a-z0-9]{4}-[a-z0-9]{4}-[a-z0-9]{12}"
UUID_REGEX = re.compile(UUID_PATTERN)


def populate_static_project(project_dir: Path, file_pattern: str = '*.json'):
    bundle_dir = project_dir / 'bundle'
    hca_dir = project_dir / 'hca'

    if hca_dir.is_dir():
        if not bundle_dir.is_dir():
            os.mkdir(str(bundle_dir))

        print(f'Removing old {file_pattern} files from: {bundle_dir}')
        for filename in glob.glob(str(bundle_dir / file_pattern)):
            new_filename = str(bundle_dir / os.path.basename(filename))
            os.unlink(new_filename)

        print(f'Populating new {file_pattern} files into: {bundle_dir}')
        for filename in glob.glob(str(hca_dir / file_pattern)):
            new_filename = str(bundle_dir / os.path.basename(filename))
            os.link(filename, new_filename)


def populate_all_static_projects(file_pattern='*.json'):
    for root, dirs, files in os.walk('projects'):
        for dir_name in dirs:
            if UUID_REGEX.match(dir_name):
                project_dir = Path('projects') / dir_name
                populate_static_project(project_dir=project_dir, file_pattern=file_pattern)
