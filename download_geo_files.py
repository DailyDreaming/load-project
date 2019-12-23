import argparse
import logging
import os
from pathlib import Path
import re
import sys
import time
from typing import (
    MutableMapping,
    Sequence
)
import uuid

import furl
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
        accession_ids = args.accessions.split(',')
        for accession_id in accession_ids:
            download_supplementary_files(accession_id)


def download_from_all_accessions():
    for sub_dir in 'existing', 'new':
        src_dir = Path('spreadsheets') / sub_dir
        paths = get_spreadsheet_paths(src_dir)
        for accession_id in [p.name.split('.')[0] for p in paths]:
            download_supplementary_files(accession_id)


def download_supplementary_files(accession_id):
    """
    Scrape web page for given accession id and download all supplementary files
    """
    logging.info('---')
    logging.info('Accession: %s', accession_id)
    project_uuid = generate_project_uuid([accession_id])

    save_file_path = f'projects/{project_uuid}/geo'
    if not os.path.exists(save_file_path):
        os.makedirs(save_file_path, exist_ok=True)
    file_link_path = f'projects/{accession_id}'
    if os.path.lexists(file_link_path):
        os.remove(file_link_path)
    os.symlink(str(project_uuid), file_link_path)

    page = requests.get(source_url_template + accession_id)

    links = supplementary_file_download_links(page.text)

    if not links:
        logging.warning('No supplementary files found on page')
        os.makedirs(save_file_path, exist_ok=True)
        open(f'{save_file_path}/no-supplementary-files', 'a').close()
        return None

    for filename, url in links:
        save_file_name = filename
        save_file_full = f'{save_file_path}/{save_file_name}'

        # download and save the file
        os.makedirs(save_file_path, exist_ok=True)
        if os.path.isfile(save_file_full):
            logging.info('Skipping existing file: %s', save_file_full)
        else:
            logging.info('Downloading to: %s', save_file_full)
            download_file(url, save_file_full)
    return save_file_path


def supplementary_file_download_links(html: str) -> Sequence:
    """
    Locate the Supplementary file table in the given html
    and pluck out the file name & urls to download the file
    """
    source = furl.furl(source_url_template)
    match = re.search(r'<table(?:(?!table).)+?>Supplementary file<(?:(?!table).)+</table>', html, re.DOTALL)
    if match:
        # For each row in the table get the filename and the http link url
        table_html = match.group()
        match = re.findall(r'<tr[^>]*>'                             # tr up to first td
                           r'<td[^>]*>([^<]+)</td>'                 # first td and capture contents
                           r'(?:(?!/tr).)*<td[^>]*>(?:(?!</td).)*'  # skip ahead another td in same tr
                           r'<a href="([^"]+)">\(http\)</a>',       # <a> tag and capture href value
                           table_html, re.DOTALL)
        if match:
            return [(filename, source.asdict()['origin'] + url if url.startswith('/') else url)
                    for filename, url in match]
        else:
            return None
    else:
        return None


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
