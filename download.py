from itertools import dropwhile
import logging
import re
import tempfile
from typing import (
    MutableMapping,
    Sequence,
    Tuple,
)

from bs4 import BeautifulSoup
import furl
from more_itertools import one
import requests

from _pathlib import Path
from create_project import (
    generate_project_uuid,
)
from util import get_target_spreadsheets

logging.basicConfig(level=logging.INFO)

source_url_template = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='


def main():
    accessions = [p.name.split('.')[0] for p in get_target_spreadsheets()]
    for accession in accessions:
        download_supplementary_files(accession)


def download_supplementary_files(accession):
    """
    Scrape web page for given accession id and download all supplementary files
    """
    logging.info('---')
    project_uuid = generate_project_uuid([accession])
    logging.info('Checking project with accession %s and UUID %s for files to download ...', accession, project_uuid)

    projects_path = Path('projects')
    geo_path = projects_path / project_uuid / 'geo'
    if not geo_path.exists():
        geo_path.mkdir(parents=True)
    create_or_update_symlink(projects_path / accession, Path(project_uuid))

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
                download_file(url, file_path)
    else:
        logging.info('No supplementary files found on %s', source_url)


def create_or_update_symlink(symlink: Path, target: Path):
    if symlink.is_symlink():
        # noinspection PyUnresolvedReferences
        current_target = symlink.readlink()
        if current_target == target:
            return
        logging.warning('Removing stale symlink from %s to %s.', symlink, current_target)
        symlink.unlink()
    elif symlink.exists():
        raise RuntimeError(f'Will not overwrite {symlink} with link to {target}')
    logging.info('Linking %s to %s', symlink, target)
    symlink.symlink_to(target, target_is_directory=True)


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


def download_file(url: str, path: Path):
    """
    Stream download the file from url, save it to path, and return response headers
    """
    with requests.get(url, stream=True) as request:
        request.raise_for_status()
        with tempfile.NamedTemporaryFile(dir=str(path.parent), delete=False) as f:
            try:
                for chunk in request.iter_content(chunk_size=1024 * 1024):
                    f.write(chunk)
            except:
                Path(f.name).unlink()
                raise
            else:
                Path(f.name).rename(path)


if __name__ == '__main__':
    main()
