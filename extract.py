from pathlib import Path
import shutil
from tarfile import TarFile
from typing import Union

import logging
logging.basicConfig(level=logging.INFO)


def main():
    # Find the RAW tar files and recursively extract them
    projects_dir = Path.cwd() / 'projects'
    for project_dir in projects_dir.iterdir():
        if project_dir.is_dir() and not project_dir.is_symlink():
            geo_dir = project_dir / 'geo'
            for geo_file in geo_dir.iterdir():
                if geo_file.name.endswith('_RAW.tar'):
                    extract_tar_file(geo_file)


def extract_tar_file(tar_path: Path):
    """
    Recursively extract tar files into a folder located at the same path as the tar file

    :param tar_path: Path to a tar file or folder containing one or more tar files
    """
    logging.debug('Running extract_tar_file(%s)', tar_path)
    if tar_path.is_dir():
        # Run this function recursively for each item found in the given dir
        logging.debug('Decending into directory %s', tar_path)
        for file in tar_path.iterdir():
            extract_tar_file(file)
    elif tar_path.is_file() and tar_path.name.endswith(('.tar', '.tar.bz', 'tar.gz')):
        # Extract the tar file and place a completion file in the new folder when done
        logging.debug('Extracting file %s', tar_path)
        file_name_stripped = file_name_minus_extensions(tar_path.name)
        if file_name_stripped and len(file_name_stripped) < len(tar_path.name):
            dest_dir = tar_path.parent / file_name_stripped
            completion_file = dest_dir / '.complete'
            if completion_file.exists():
                logging.info('Expansion of %s already complete', tar_path)
            else:
                with TarFile.open(tar_path.as_posix()) as tar_file:
                    assert completion_file.name not in tar_file.getnames()
                    if dest_dir.exists():
                        logging.info('Removing partially expanded %s', dest_dir)
                        shutil.rmtree(dest_dir.as_posix())
                    logging.info('Expanding %s', dest_dir)
                    dest_dir.mkdir()
                    tar_file.extractall(dest_dir.as_posix())
                    completion_file.touch()
                    logging.info('Expansion of %s is complete', dest_dir)
            # Run this function on the new folder where the tar file was extracted
            extract_tar_file(dest_dir)


def file_name_minus_extensions(file_name: str) -> Union[str, None]:
    """
    Return the file name minus the extension(s)

    >>> file_name_minus_extensions('foo.tar')
    foo
    >>> file_name_minus_extensions('foo.tar.gz')
    foo
    """
    if '.' in file_name:
        return file_name[:file_name.index('.')]
    else:
        return None


if __name__ == '__main__':
    main()
