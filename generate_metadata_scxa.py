from typing import List

from _pathlib import Path

from more_itertools import one

from create_project import timestamp, write_project_json
from util import get_target_project_dirs, generate_project_uuid


def main():
    project_dirs = get_target_project_dirs()

    for project_dir in project_dirs:
        generate_metadata(project_dir)


def generate_metadata(project_dir: Path):
    project_json = create_project_json(project_dir)
    write_project_json(project_json, bundle_dir(project_dir))


def bundle_dir(project_dir: Path):
    return project_dir / 'bundle'


def create_project_json(project_dir: Path):
    accession = accession_from_project_dir(project_dir)
    idf_file = file_by_suffix(scxa_dir(project_dir), 'idf.txt')
    idf_metadata = parse_mage_tab(idf_file)
    version = timestamp()
    return {
        "describedBy": "https://schema.dev.data.humancellatlas.org/type/project/14.0.0/project",
        "schema_type": "project",
        "project_core": {
            "project_short_name": only(idf_metadata['Investigation Title']),
            "project_title": only(idf_metadata['Investigation Title']),
            "project_description": only(idf_metadata['Experiment Description']),
        },
        "provenance": {
            "document_id": generate_project_uuid(accession),
            "submission_date": version,
            "update_date": version
        }
    }


def accession_from_project_dir(project_dir):
    accession = project_dir.name
    assert accession.startswith('E-')
    return accession


def only(row: List[str]):
    """
    Return only the first element and ensure any others are ''
    >>> only(['foo', '', '', ''])
    'foo'
    >>> only(['foo', '', '', 'bar'])
    Traceback (most recent call last):
        ...
    AssertionError
    """
    first, *rest = row
    assert all(cell == '' for cell in rest)
    return first


def parse_mage_tab(idf_file: Path):
    data = {}
    with idf_file.open('r') as f:
        assert f.readline().startswith('MAGE-TAB Version')
        for line in f.readlines():
            line = line.rstrip('\n')
            cols = line.split('\t')
            if len(cols) > 0:
                header, *values = cols
                data[header] = values
    return data


def file_by_suffix(scxa: Path, suffix: str):
    metadata_dir = scxa / 'experiment-metadata'
    return one(file for file in metadata_dir.iterdir() if file.name.endswith(suffix))


def scxa_dir(project_dir):
    return project_dir / 'scxa'
