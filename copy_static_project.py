import os
import glob
from pathlib import Path


def populate_static_project(project_dir: Path, file_pattern: str = '*.json'):
    assert project_dir.is_dir()
    assert len(project_dir.name) == len('4a95101c-9ffc-4f30-a809-f04518a23803')
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
    else:
        raise RuntimeError(f'{hca_dir} does not exist.  Nothing to link over.')


def populate_all_static_projects(file_pattern: str = '*.json'):
    for root, dirs, files in os.walk('projects'):
        for project_uuid in dirs:
            if len(project_uuid) == len('4a95101c-9ffc-4f30-a809-f04518a23803'):
                src_dir = Path('projects') / project_uuid
                hca_dir = src_dir / 'hca'
                if hca_dir.is_dir():
                    populate_static_project(src_dir, file_pattern=file_pattern)
                    print(f'Hard-linked project ("{file_pattern}" only) contents: '
                          f'{project_uuid}/hca to {project_uuid}/bundle.')
