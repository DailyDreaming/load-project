import argparse
from dataclasses import dataclass
import logging
from _pathlib import Path
import sys
from typing import (
    Mapping,
    Sequence,
)
from uuid import UUID

from count_cells import CountCells
from create_project import (
    generate_project_uuid,
)
from util import (
    get_target_project_dirs,
    get_target_spreadsheets,
)

logging.basicConfig(level=logging.INFO)

projects_path = Path('projects')


@dataclass
class ProjectReport:
    uuid: UUID = None
    accession: str = None
    project_path: Path = None  # projects/{uuid}
    symlink: Path = None  # projects/{accession} symlink to project_path
    spreadsheet: Path = None  # spreadsheets/(new|existing)/{accession}.0.xlsx
    geo_files: int = 0  # number of downloaded geo files in projects/{uuid}/geo
    num_matrices: int = 0  # number of matrices in projects/{uuid}/matrices
    zipped_matrix: Path = None  # projects/{uuid}/bundle/matrix.mtx.zip
    cell_count: int = 0  # number of cells counted
    gene_count: int = 0  # number of genes counted
    num_metadata_files: int = 0  # number of metadata JSON files in projects/{uuid}/bundle
    num_hca_metadata_files: int = 0  # number of metadata JSON files in projects/{uuid}/hca
    zipped_hca_matrix = None  # projects/{uuid}/hca/matrix.mtx.zip

    @classmethod
    def tsv_header(cls) -> Sequence[str]:
        return [
            'uuid',
            'accession',
            'project_path',
            'symlink',
            'spreadsheet',
            'geo_files',
            'num_matrices',
            'zipped_matrix',
            'cell_count',
            'gene_count',
            'num_metadata_files',
            'num_hca_metadata_files',
            'zipped_hca_matrix',
        ]

    def tsv_row(self) -> Sequence[str]:
        return ['-' if x is None else str(x) for x in [
            self.uuid,
            self.accession,
            self.project_path,
            self.symlink,
            self.spreadsheet,
            self.geo_files,
            self.num_matrices,
            self.zipped_matrix,
            self.cell_count,
            self.gene_count,
            self.num_metadata_files,
            self.num_hca_metadata_files,
            self.zipped_hca_matrix
        ]]


def write_report_as_tsv(report: Mapping[UUID, ProjectReport]):
    """
    Write a report as a tsv to stdout

    :param report: A dict mapping UUIDs to objects containing project details
    """
    logging.debug('Starting tsv output')
    print('\t'.join(ProjectReport.tsv_header()))
    for uuid, project in report.items():
        print('\t'.join(project.tsv_row()))
    logging.debug('Finished tsv output')


def overview_report() -> Mapping[UUID, ProjectReport]:
    """
    Generate a report that reconciles the presence of certain resources associated with each project

    :return: An overview report in the form of a dict mapping UUIDs to objects containing project details
    """

    report = {}

    logging.debug('Searching for project uuids in the projects path ...')
    for uuid in get_project_uuids():
        project_path = projects_path / str(uuid)
        report[uuid] = ProjectReport(uuid=uuid,
                                     project_path=project_path)

    logging.debug('Searching for accession ids in the projects path ...')
    for accession_id in [p.name for p in get_target_project_dirs()]:
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
                report[uuid_from_symlink].accession = accession_id
                report[uuid_from_symlink].symlink = accession_symlink
            else:  # Create a new project in the report using the uuid from symlink's target
                assert not (projects_path / str(uuid_from_symlink)).exists(), \
                    f'get_project_uuids() failed to find {str(uuid_from_symlink)}'
                logging.debug('New accession %s found as accession with symlink, adding uuid %s',
                              accession_id, str(uuid_from_symlink))
                report[uuid_from_symlink] = ProjectReport(uuid=uuid_from_symlink,
                                                          accession=accession_id,
                                                          symlink=accession_symlink)
        else:
            uuid_from_accession_id = UUID(generate_project_uuid(accession_id))
            if uuid_from_accession_id in report:  # update existing project in report
                report[uuid_from_accession_id].accession = accession_id
            else:
                logging.debug('New accession %s found as accession without symlink, adding uuid %s',
                              accession_id, str(uuid_from_accession_id))
                report[uuid_from_accession_id] = ProjectReport(uuid=uuid_from_accession_id,
                                                               accession=accession_id)

    logging.debug('Searching for spreadsheets ...')
    for accession_id, file in get_target_spreadsheets():
        logging.debug('Checking: %s', file)
        uuid = UUID(generate_project_uuid(accession_id))
        try:
            report[uuid].spreadsheet = file
        except KeyError:
            logging.debug('New accession %s found in spreadsheets, adding uuid %s',
                          accession_id, str(uuid))
            report[uuid] = ProjectReport(uuid=uuid,
                                         accession=accession_id,
                                         spreadsheet=file)

    logging.debug('Fetching cell count ...')
    for accession_id, cell_count in CountCells.get_cached_cell_counts().items():
        uuid = UUID(generate_project_uuid(accession_id))
        try:
            report[uuid].cell_count = cell_count
        except KeyError:
            logging.debug('New accession %s found, adding uuid %s',
                          accession_id, str(uuid))
            report[uuid] = ProjectReport(uuid=uuid,
                                         accession=accession_id,
                                         cell_count=cell_count)

    logging.debug('Fetching gene count ...')
    for accession_id, gene_count in CountCells.get_cached_gene_counts().items():
        uuid = UUID(generate_project_uuid(accession_id))
        try:
            report[uuid].gene_count = gene_count
        except KeyError:
            logging.debug('New accession %s found, adding uuid %s',
                          accession_id, str(uuid))
            report[uuid] = ProjectReport(uuid=uuid,
                                         accession=accession_id,
                                         gene_count=gene_count)

    logging.debug('Counting geo files ...')
    for uuid in report:
        path = projects_path / str(uuid) / 'geo'
        report[uuid].geo_files = get_file_count(path, glob='**/*')

    logging.debug('Counting matrices ...')
    for uuid in report:
        path = projects_path / str(uuid) / 'matrices'
        report[uuid].num_matrices = get_file_count(path, glob='**/matrix.mtx.gz')

    logging.debug('Checking for zipped_matrix ...')
    for uuid in report:
        zipped_matrix = projects_path / str(uuid) / 'bundle' / 'matrix.mtx.zip'
        if zipped_matrix.is_file():
            report[uuid].zipped_matrix = zipped_matrix

    logging.debug('Checking for metadata_json_count ...')
    for uuid in report:
        path = projects_path / str(uuid) / 'bundle'
        report[uuid].num_metadata_files = get_file_count(path, glob='*.json')

    logging.debug('Checking for num_hca_metadata_files ...')
    for uuid in report:
        path = projects_path / str(uuid) / 'hca'
        report[uuid].num_hca_metadata_files = get_file_count(path, glob='*.json')

    logging.debug('Checking for zipped_matrix ...')
    for uuid in report:
        zipped_matrix = projects_path / str(uuid) / 'hca' / 'matrix.mtx.zip'
        if zipped_matrix.is_file():
            report[uuid].zipped_hca_matrix = zipped_matrix

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
    write_report_as_tsv(overview_report())
    # print("Done.")


if __name__ == '__main__':
    main(sys.argv[1:])
