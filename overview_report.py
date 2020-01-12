import argparse
import csv
from dataclasses import dataclass
import logging
from pathlib import Path
import sys
from typing import (
    Mapping,
    Sequence,
)
from uuid import UUID

from count_cells import (
    get_accession_ids,
    get_cell_counts,
)
from create_project import (
    generate_project_uuid,
    get_spreadsheet_paths,
)

logging.basicConfig(level=logging.INFO)

projects_path = Path('projects')


@dataclass
class ProjectReport:
    uuid: UUID = None
    accession_id: str = None
    project_path: Path = None  # projects/{uuid}
    accession_symlink: Path = None  # projects/{accession} symlink to project_path
    spreadsheet: Path = None  # spreadsheets/(new|existing)/{accession}.0.xlsx
    geo_files_count: int = 0  # number of downloaded geo files in projects/{uuid}/geo
    matrices_count: int = 0  # number of matrices in projects/{uuid}/matrices
    zipped_matrix: Path = None  # projects/{uuid}/bundle/matrix.mtx.zip
    cell_count: int = 0  # value of cell count in cell_counts.json
    metadata_json_count: int = 0  # number of metadata JSON files in projects/{uuid}/bundle

    @classmethod
    def tsv_header(cls) -> Sequence[str]:
        return [
            'UUID',
            'Accession ID',
            'Project path',
            'Accession symlink',
            'Spreadsheet',
            'Geo files count',
            'Matrices count',
            'Zipped matrix',
            'Cell count',
            'Metadata JSON count',
        ]

    def tsv_row(self) -> Sequence[str]:
        return [
            self.uuid,
            self.accession_id,
            self.project_path,
            self.accession_symlink,
            self.spreadsheet,
            self.geo_files_count,
            self.matrices_count,
            self.zipped_matrix,
            self.cell_count,
            self.metadata_json_count,
        ]


def write_report_to_tsv(report: Mapping[UUID, ProjectReport], output_file: Path):
    """
    Write a report to a tsv file

    :param report: A dict mapping UUIDs to objects containing project details
    :param output_file: A path to the desired output file
    """
    logging.debug('Opening %s for tsv writing', str(output_file))
    with open(str(output_file), 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        logging.debug('Writting tsv header row')
        writer.writerow(ProjectReport.tsv_header())
        for uuid, project in report.items():
            logging.debug('Writting tsv row for uuid %s', str(uuid))
            writer.writerow(project.tsv_row())
    logging.debug('Closed %s', str(output_file))


def overview_report() -> Mapping[UUID, ProjectReport]:
    """
    Generate a report that reconciles the presence of certain resources associated with each project

    :return: An overview report in the form of a dict mapping UUIDs to objects containing project details
    """

    report = {}

    # ---
    logging.debug('Searching for project uuids in the projects path ...')
    for uuid_from_symlink in get_project_uuids():
        project_path = projects_path / str(uuid_from_symlink)
        report[uuid_from_symlink] = ProjectReport(uuid=uuid_from_symlink,
                                                  project_path=project_path)

    # ---
    logging.debug('Searching for accession ids in the projects path ...')
    for accession_id in get_accession_ids():
        # accession ids are symlinks to folders named by uuid
        accession_symlink = projects_path / accession_id
        assert accession_symlink.is_symlink()
        expanded_path = accession_symlink.resolve()
        if not expanded_path.exists():
            logging.debug('Error: Symlink %s has invalid target %s', accession_id, expanded_path)
            expanded_path = None
        try:
            uuid_from_symlink = UUID(expanded_path.name)
        except AttributeError:
            # logging.debug('Error: UUID(None.name)')
            uuid_from_symlink = None  # Symlink has an invalid target
        except ValueError:
            logging.debug('Error: Symlink %s target %s is invalid UUID', accession_id, expanded_path.name)
            uuid_from_symlink = None  # Value from symlink target wasn't a valid UUID
        if uuid_from_symlink:
            if uuid_from_symlink in report:  # Update existing project in report with accession info
                report[uuid_from_symlink].accession_id = accession_id
                report[uuid_from_symlink].accession_symlink = accession_symlink
            else:  # Create a new project in the report using the uuid from symlink's target
                assert not (projects_path / str(uuid_from_symlink)).exists(), \
                    f'get_project_uuids() failed to find {str(uuid_from_symlink)}'
                logging.debug('New accession %s found as accession with symlink, adding uuid %s',
                              accession_id, str(uuid_from_symlink))
                report[uuid_from_symlink] = ProjectReport(uuid=uuid_from_symlink,
                                                          accession_id=accession_id,
                                                          accession_symlink=accession_symlink)
        else:
            uuid_from_accession_id = UUID(generate_project_uuid(accession_id))
            if uuid_from_accession_id in report:  # update existing project in report
                report[uuid_from_accession_id].accession_id = accession_id
            else:
                logging.debug('New accession %s found as accession without symlink, adding uuid %s',
                              accession_id, str(uuid_from_accession_id))
                report[uuid_from_accession_id] = ProjectReport(uuid=uuid_from_accession_id,
                                                               accession_id=accession_id)

    # ---
    logging.debug('Searching for spreadsheets ...')
    for sub_dir in 'existing', 'new':
        for file in get_spreadsheet_paths(Path(f'spreadsheets/{sub_dir}')):
            logging.debug('Checking: %s', file)
            accession_id = file.name[:-len('.0.xlsx')]
            uuid_from_symlink = UUID(generate_project_uuid(accession_id))
            try:
                report[uuid_from_symlink].spreadsheet = file
            except KeyError:
                logging.debug('New accession %s found in spreadsheets, adding uuid %s',
                              accession_id, str(uuid_from_symlink))
                report[uuid_from_symlink] = ProjectReport(uuid=uuid_from_symlink,
                                                          accession_id=accession_id,
                                                          spreadsheet=file)

    # ---
    logging.debug('Reading cell_counts.json ...')
    for accession_id, cell_count in get_cell_counts().items():
        uuid_from_symlink = UUID(generate_project_uuid(accession_id))
        try:
            report[uuid_from_symlink].cell_count = cell_count
        except KeyError:
            logging.debug('New accession %s found in cell_counts.json, adding uuid %s',
                          accession_id, str(uuid_from_symlink))
            report[uuid_from_symlink] = ProjectReport(uuid=uuid_from_symlink,
                                                      accession_id=accession_id,
                                                      cell_count=cell_count)

    # ---
    logging.debug('Counting geo files ...')
    for uuid_from_symlink in report:
        path = projects_path / str(uuid_from_symlink) / 'geo'
        report[uuid_from_symlink].geo_files_count = get_file_count(path, glob='**/*')

    # ---
    logging.debug('Counting matrices ...')
    for uuid_from_symlink in report:
        path = projects_path / str(uuid_from_symlink) / 'matrices'
        report[uuid_from_symlink].matrices_count = get_file_count(path, glob='**/matrix.mtx.gz')

    # ---
    logging.debug('Checking for zipped_matrix ...')
    for uuid_from_symlink in report:
        zipped_matrix = projects_path / str(uuid_from_symlink) / 'bundle' / 'matrix.mtx.zip'
        if zipped_matrix.is_file():
            report[uuid_from_symlink].zipped_matrix = zipped_matrix

    # ---
    logging.debug('Checking for metadata_json_count ...')
    for uuid_from_symlink in report:
        path = projects_path / str(uuid_from_symlink) / 'bundle'
        report[uuid_from_symlink].metadata_json_count = get_file_count(path, glob='*.json')

    return report


def get_file_count(path: Path, glob: str) -> int:
    """
    Return the count of files in a folder

    :param path: A path to a folder to check
    :param glob: The glob pattern

    :return: Number of files counted
    """
    if path.exists() and path.is_dir():
        return sum([1 for f in path.glob(glob) if f.is_file()])
    else:
        return 0


def get_project_uuids() -> Sequence[UUID]:
    """
    Return a list of all the project UUIDs found in the projects path
    """
    uuids = []
    for p in projects_path.iterdir():
        logging.debug('Checking: %s', p)
        if p.is_dir() and not p.is_symlink():
            try:
                uuid = UUID(p.name)
            except ValueError:
                continue  # skip value if not a valid UUID
            logging.debug('Found: %s', uuid)
            uuids.append(uuid)
    return uuids


def main(argv):
    """
    Support for command line execution
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verbose', '-v',
                        action='store_true',
                        help='Verbose debug output')
    args = parser.parse_args(argv)
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    write_report_to_tsv(overview_report(), Path('overview_report.tsv'))
    print("File written: ./overview_report.tsv")


if __name__ == '__main__':
    main(sys.argv[1:])
