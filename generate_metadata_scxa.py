import json
import logging
from typing import List

from more_itertools import one

from _pathlib import Path
from count_cells import CountCells
from create_project import (
    generate_analysis_json,
    generate_links_json,
    timestamp,
    write_project_json,
)
from util import (
    generate_file_uuid,
    generate_project_uuid,
    get_target_project_dirs,
)

log = logging.getLogger(__name__)


def main():
    for project_dir in get_target_project_dirs():
        accession = project_dir.name
        if accession.startswith('E-'):
            generate_metadata(accession, project_dir)


def generate_metadata(accession: str, project_dir: Path):
    project_uuid = generate_project_uuid(accession)
    assert project_dir.readlink().name == project_uuid
    bundle_dir = project_dir / 'bundle'
    project_json = create_project_json(accession, project_dir, project_uuid)
    write_project_json(project_json, bundle_dir)
    generate_analysis_json(project_uuid, bundle_dir)
    cell_count = CountCells.get_cached_cell_count(project_dir)
    generate_cell_suspension_json(bundle_dir, cell_count, project_uuid)
    generate_links_json(bundle_dir)


def generate_cell_suspension_json(bundle_dir, cell_count, bundle_uuid):
    file_name = 'cell_suspension_0.json'
    file_uuid = generate_file_uuid(bundle_uuid, file_name)
    cell_json = create_cell_suspension_json(cell_count, file_uuid)
    output_file = bundle_dir / file_name
    with open(output_file, 'w') as f:
        json.dump(cell_json, f, indent=4)
    print(f'"{output_file}" successfully written.')


def create_cell_suspension_json(cell_count, file_uuid):
    version = timestamp()
    cell_suspension_json = {
        "describedBy": "https://schema.dev.data.humancellatlas.org/type/biomaterial/13.1.0/cell_suspension",
        "schema_type": "biomaterial",
        "biomaterial_core": {
            "biomaterial_id": str(file_uuid),
            "biomaterial_name": None,
            "biomaterial_description": None,
            "ncbi_taxon_id": [0],
            "genotype": None
        }, "genus_species": [
            {
                "text": None,
                "ontology": None,
                "ontology_label": None
            }
        ], "provenance": {
            "document_id": str(file_uuid),
            "submission_date": version,
            "update_date": version
        },
        "estimated_cell_count": cell_count,
    }
    return cell_suspension_json


def create_project_json(accession: str, project_dir: Path, project_uuid):
    idf_file = file_by_suffix(scxa_dir(project_dir), 'idf.txt')
    idf_metadata = parse_mage_tab(idf_file)
    version = timestamp()
    return {
        "describedBy": "https://schema.dev.data.humancellatlas.org/type/project/14.0.0/project",
        "schema_type": "project",
        "project_core": {
            "project_short_name": accession,
            "project_title": only(idf_metadata['Investigation Title']),
            **({
                   "project_description": only(idf_metadata['Experiment Description']),
               } if 'Experiment Description' in idf_metadata else {})
        },
        "provenance": {
            "document_id": project_uuid,
            "submission_date": version,
            "update_date": version
        }
    }


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
    log.debug('Parsing IDF file "%s"', idf_file)
    with idf_file.open('r') as f:
        for line in f.readlines():
            line = line.rstrip('\n')
            cols = line.split('\t')
            if len(cols) > 0:
                header, *values = cols
                data[header] = values
    assert 'MAGE-TAB Version' in data
    return data


def file_by_suffix(scxa: Path, suffix: str):
    metadata_dir = scxa / 'experiment-metadata'
    return one(file for file in metadata_dir.iterdir() if file.name.endswith(suffix))


def scxa_dir(project_dir):
    return project_dir / 'scxa'


if __name__ == "__main__":
    logging.basicConfig(format='%(asctime)s %(levelname)s %(threadName)s: %(message)s', level=logging.DEBUG)
    main()
