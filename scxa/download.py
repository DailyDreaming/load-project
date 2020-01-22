import logging
import os
import sys
import uuid
from collections import defaultdict
from concurrent.futures.thread import ThreadPoolExecutor
from pathlib import Path
from typing import NamedTuple, List
from urllib.request import urlretrieve

import requests

log = logging.getLogger(__name__)


namespace_uuid = '1edd6b61-fd3d-4dee-b911-958413871a5c'


def download_projects(path: Path):
    accessions = all_accessions()
    for accession in accessions:
        download_project(accession, path)


def all_accessions():
    projects = list_projects()
    accessions = [p['experimentAccession'] for p in projects]
    return accessions


def download_projects_parallel(path: Path):
    with ThreadPoolExecutor() as executor:
        for accession in all_accessions():
            accessions_to_futures = defaultdict(set)
            scxa_path = make_and_link_download_dir(accession, path)
            files = project_files(accession)
            for file in files:
                accessions_to_futures[accession].add(executor.submit(file.idempotent_download, scxa_path))
    failed_projects = []
    for accession, futures in accessions_to_futures.items():
        if not all(future.result() for future in futures):
            failed_projects.append(accession)

    print('Failed projects', file=sys.stderr)
    for project in failed_projects:
        print(project, file=sys.stderr)


def download_project(accession: str, path: Path):
    scxa_path = make_and_link_download_dir(accession, path)
    files = project_files(accession)
    for file in files:
        file.idempotent_download(scxa_path)


def project_files(accession):
    return [
        FileURL(accession=accession, file_type=file_type, zipped=zipped)
        for file_type, zipped in [
            ('experiment-metadata', True),
            ('experiment-design', False),
            ('cluster', False),
            ('marker-genes', True),
            ('normalised', True),
            ('quantification-raw', True),
        ]
    ]


def make_and_link_download_dir(accession, path):
    uuid_ = uuid.uuid5(uuid.UUID(namespace_uuid), accession)
    project_dir = path / str(uuid_)
    project_dir.mkdir(exist_ok=True)
    try:
        os.symlink(project_dir.absolute(), path / accession)
    except FileExistsError:
        pass
    scxa_path = project_dir / 'scxa'
    scxa_path.mkdir(exist_ok=True)
    return scxa_path


def list_projects() -> List[dict]:
    projects_url = 'https://www.ebi.ac.uk/gxa/sc/json/experiments'
    response = requests.get(projects_url)
    return response.json()['experiments']


class FileURL(NamedTuple):
    accession: str
    file_type: str
    zipped: bool

    def to_url(self):
        base_url = 'https://www.ebi.ac.uk/gxa/sc/experiment'
        return f"{base_url}/{self.accession}/download{'/zip' if self.zipped else ''}?fileType={self.file_type}"

    def idempotent_download(self, path) -> bool:
        name = self.file_type + ('.zip' if self.zipped else '')
        file_path = path / name
        if not file_path.exists():
            log.info('Downloading new file `%s` from URL `%s`', file_path, self.to_url())
            try:
                urlretrieve(self.to_url(), filename=file_path)
                return True
            except Exception:
                log.warning('Failed to download file `%s` from URL `%s`', exc_info=True)
                return False
        else:
            log.info('Skipping download of file `%s`', file_path)
            return True


if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s', level=logging.INFO)

    download_projects_parallel(Path('projects'))
    # download_projects(Path('test'))
