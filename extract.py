from operator import methodcaller
from _pathlib import Path
import shutil
from tarfile import TarFile

import logging

from util import get_target_project_dirs

logging.basicConfig(level=logging.INFO)


def main():

    project_dirs = get_target_project_dirs()

    for project_dir in project_dirs:
        extract_tar_file_recursive(project_dir)


def extract_tar_file(tar_path: Path, dest_path: Path):
    """
    Extract a tar file  and put completion file in destination folder once complete.
    Skips extraction if a completion file found in the destination folder

    :param tar_path: Path to a tar file
    :param dest_path: Path to put extracted tar contents
    """
    completion_file = dest_path / '.complete'
    if completion_file.exists():
        logging.info('Expansion of %s already complete', tar_path)
    else:
        openmode = 'r:gz' if tar_path.name.endswith('.tar.gz') else 'r'
        with TarFile.open(str(tar_path), mode=openmode) as tar_file:
            assert completion_file.name not in tar_file.getnames()
            if dest_path.exists():
                logging.info('Removing partially expanded %s', dest_path)
                shutil.rmtree(str(dest_path))
            logging.info('Expanding %s', dest_path)
            dest_path.mkdir()
            tar_file.extractall(str(dest_path))
            completion_file.touch()
            logging.info('Expansion of %s is complete', dest_path)


def extract_tar_file_recursive(tar_path: Path):
    """
    Recursively extract tar files into a folder located at the same path as the tar file

    :param tar_path: Path to a tar file or folder containing one or more tar files
    """
    logging.debug('Running extract_tar_file_recursive(%s)', tar_path)
    if tar_path.is_dir():
        logging.debug('Decending into directory %s', tar_path)
        # Iterate over directory contents sorted by type with files first, then
        # directories. This is done to avoid wasting time processing a directory
        # that could be itself be deleted and re-extracted from a tar file
        for file in sorted(tar_path.iterdir(), key=methodcaller('is_dir')):
            extract_tar_file_recursive(file)
    elif tar_path.name.endswith(('.tar', '.tar.gz')):
        assert tar_path.is_file()
        logging.debug('Extracting file %s', tar_path)
        base_file_name = get_base_file_name(tar_path.name)
        assert 0 < len(base_file_name) < len(tar_path.name)
        dest_path = tar_path.parent / base_file_name
        extract_tar_file(tar_path, dest_path)  # Extract the tar file to a subfolder
        extract_tar_file_recursive(dest_path)  # Check subfolder for tar files and extract them


def get_base_file_name(file_name: str) -> str:
    """
    Return the file name minus the extension(s)

    >>> get_base_file_name('foo.tar')
    foo
    >>> get_base_file_name('foo.tar.gz')
    foo
    """
    pos = file_name.index('.')
    if pos:
        return file_name[:pos]
    else:
        return file_name


if __name__ == '__main__':
    main()
