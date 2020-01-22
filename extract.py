import zipfile
from operator import methodcaller
from _pathlib import Path
import shutil
import tarfile

import logging

from util import get_target_project_dirs

logging.basicConfig(level=logging.INFO)


def main():
    project_dirs = get_target_project_dirs()

    for project_dir in project_dirs:
        extract_recursive(project_dir)


def extract_file(src_path: Path, dest_path: Path, compression='tar'):
    """
    Extract a compressed file and put completion file in destination folder once complete.
    Skips extraction if a completion file found in the destination folder

    :param src_path: Path to a compressed file
    :param dest_path: Path to put extracted contents
    :param compression: Either 'tar' or 'zip'
    """
    completion_file = dest_path / '.complete'
    if completion_file.exists():
        logging.info('Expansion of %s already complete', src_path)
    else:
        if compression == 'tar':
            openmode = 'r:gz' if src_path.name.endswith('.tar.gz') else 'r'
            extractor = tarfile.TarFile.open(str(src_path), mode=openmode)
            assert completion_file.name not in extractor.getnames()
        elif compression == 'zip':
            extractor = zipfile.ZipFile(str(src_path), 'r')
        else:
            raise ValueError('Unsupported compression')
        with extractor:
            if dest_path.exists():
                logging.info('Removing partially expanded %s', dest_path)
                shutil.rmtree(str(dest_path))
            logging.info('Expanding %s', dest_path)
            dest_path.mkdir()
            extractor.extractall(str(dest_path))
            completion_file.touch()
            logging.info('Expansion of %s is complete', dest_path)


def extract_recursive(compressed_path: Path):
    """
    Recursively extract tar files into a folder located at the same path as the tar file

    :param compressed_path: Path to a compressed file or folder containing one or more tar files
    """
    logging.debug('Running extract_recursive(%s)', compressed_path)
    if compressed_path.is_dir():
        logging.debug('Decending into directory %s', compressed_path)
        # Iterate over directory contents sorted by type with files first, then
        # directories. This is done to avoid wasting time processing a directory
        # that could be itself be deleted and re-extracted from a tar file
        for file in sorted(compressed_path.iterdir(), key=methodcaller('is_dir')):
            extract_recursive(file)
    elif compressed_path.name.endswith(('.tar', '.tar.gz')):
        assert compressed_path.is_file()
        logging.debug('Extracting file %s', compressed_path)
        base_file_name = get_base_file_name(compressed_path.name)
        assert 0 < len(base_file_name) < len(compressed_path.name)
        dest_path = compressed_path.parent / base_file_name
        extract_file(compressed_path, dest_path, compression='tar')  # Extract the tar file to a subfolder
        extract_recursive(dest_path)  # Check subfolder for tar files and extract them
    elif compressed_path.name in ('experiment-metadata.zip',
                                  'experiment-metadata.zip',
                                  'marker-genes.zip',
                                  'normalised.zip',
                                  'quantification-raw.zip'):
        # This is a zip download from SCXA
        base_file_name = get_base_file_name(compressed_path.name)
        dest_path = compressed_path.parent / base_file_name
        extract_file(compressed_path, dest_path, compression='zip')


def get_base_file_name(file_name: str) -> str:
    """
    Return the file name minus the extension(s)

    >>> get_base_file_name('foo.tar')
    'foo'
    >>> get_base_file_name('foo.tar.gz')
    'foo'
    >>> get_base_file_name('marker-genes.zip')
    'marker-genes'
    """
    pos = file_name.index('.')
    if pos:
        return file_name[:pos]
    else:
        return file_name


if __name__ == '__main__':
    main()
