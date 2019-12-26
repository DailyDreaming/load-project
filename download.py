import argparse
from itertools import dropwhile
import logging
import os
from pathlib import Path
import re
import sys
from typing import (
    MutableMapping,
    Sequence,
    Tuple,
)

from bs4 import BeautifulSoup
import furl
from more_itertools import one
import requests

from create_project import (
    generate_project_uuid,
    get_spreadsheet_paths,
)

logging.basicConfig(level=logging.INFO)

source_url_template = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--all', '-A',
                       action='store_true')
    group.add_argument('--accessions', '-a',
                       help='Comma separated list of GEO accession ids')
    args = parser.parse_args(argv)

    if args.all:
        download_from_all_accessions()
    else:
        accessions = args.accessions.split(',')
        for accession in accessions:
            download_supplementary_files(accession)


def download_from_all_accessions():
    for sub_dir in 'existing', 'new':
        src_dir = Path('spreadsheets') / sub_dir
        paths = get_spreadsheet_paths(src_dir)
        for accession in [p.name.split('.')[0] for p in paths]:
            download_supplementary_files(accession)


def download_supplementary_files(accession):
    """
    Scrape web page for given accession id and download all supplementary files
    """
    logging.info('---')
    project_uuid = generate_project_uuid([accession])
    logging.info('Downloading files for project accession %s, UUID %s.', accession, project_uuid)

    projects_path = Path.cwd() / 'projects'
    project_path = projects_path / project_uuid
    geo_path = project_path / 'geo'
    if not geo_path.exists():
        geo_path.mkdir(parents=True)
    accession_path = projects_path / accession
    if not accession_path.is_symlink() or accession_path.resolve() != project_path.resolve():
        logging.info('Linking %s to %s', accession_path, project_path)
        if accession_path.exists():
            accession_path.unlink()
        accession_path.symlink_to(project_path, target_is_directory=True)

    source_url = source_url_template + accession
    page = requests.get(source_url)
    links = supplementary_file_download_links(accession, page.text)
    if links:
        for file_name, url in links:
            file_path = geo_path / file_name
            if file_path.is_file():
                logging.info('Skipping existing file: %s', file_path)
            else:
                logging.info('Downloading to: %s', file_path)
                download_file(url, file_path.as_posix())
    else:
        logging.info('No supplementary files found on %s', source_url)


def supplementary_file_download_links(accession, html: str) -> Sequence[Tuple[str, str]]:
    """
    Locate the supplementary file table in the given html and pluck out the
    file name & URLs to download the file.
    """
    source = furl.furl(source_url_template)
    base_url = source.asdict()['origin']
    html = BeautifulSoup(html, 'html.parser')
    file_name_column_name = 'Supplementary file'
    url_column_name = 'Download'
    # Find single table listing supplementary files
    tables = [td.parent.parent for td in html('td', string=file_name_column_name)]
    if tables:
        table = one(tables)
        # Determine file name and URL columns
        table_header, *table_body = table('tr')
        file_name_column_index = table_header('td').index(one(table_header('td', string=file_name_column_name)))
        url_column_index = table_header('td').index(one(table_header('td', string=url_column_name)))
        # Drop single-cell rows with notes at the bottom of the table
        table_body = reversed(list(dropwhile(lambda tr: len(tr('td')) == 1, reversed(table_body))))
        links = []
        for tr in table_body:
            td = tr('td')
            file_name = td[file_name_column_index].string
            url = td[url_column_index].find('a', text='(http)').attrs['href']
            links.append((file_name, base_url + url if url.startswith('/') else url))
        return links
    else:
        # Otherwise find a table cell explcit declaring absence of files
        if html('td', string='Supplementary data files not provided'):
            return []
        else:
            # Otherwise find error message delaring accession as inaccessible
            pattern = re.compile(f'Accession "{accession}" is currently private and is scheduled to be released on')
            if html(string=pattern):
                return []
            else:
                assert False, "No known method for extracting files from GEO listing applies"


def filename_from_headers(headers: MutableMapping[str, str]):
    """
    Get the filename from the headers
    """
    if isinstance(headers, MutableMapping) and 'Content-Disposition' in headers:
        match = re.search(r'filename="([^"]+)"', headers['Content-Disposition'])
        if match:
            return match.group(1)
    return None


def download_file(url: str, path: str) -> dict:
    """
    Stream download the file from url, save it to path, and return response headers
    """
    try:
        with requests.get(url, stream=True) as request:
            request.raise_for_status()
            with open(path, 'wb') as f:
                for chunk in request.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
            headers = request.headers
        return headers
    except KeyboardInterrupt:
        logging.warning('Download canceled ...')
        os.remove(path)
        return {}


if __name__ == '__main__':
    main(sys.argv[1:])
