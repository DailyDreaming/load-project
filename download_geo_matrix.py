import argparse
import furl
import os
import re
import requests
import sys
import time
from typing import MutableMapping
import uuid

from create_project import generate_project_uuid

import logging
logging.basicConfig(level=logging.INFO)

source_url_template = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='


def main(argv):

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--accessions', '-a',
                        required=True,
                        help='Comma separated list of GEO accession ids')
    args = parser.parse_args(argv)

    accessions = args.accessions.split(',')

    for acc in accessions:
        download_supplementary_files(acc)


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

    source = furl.furl(source_url_template)

    page = requests.get(source_url_template + accession_id)

    if 'Supplementary file' not in page.text:
        logging.warning('No supplementary files found on page')
        os.makedirs(save_file_path, exist_ok=True)
        open(f'{save_file_path}/no-supplementary-files', 'a').close()
        return None

    # extract the html links
    html = page.text[page.text.rindex('Supplementary file'):]
    links = re.findall(r'<a href="[^"]+">\(http\)</a>', html)
    logging.info('Found %s links', len(links))

    for link in links:
        # extract the urls from the html links
        match = re.search(r'href="([^"]+)"', link)
        if not match:
            logging.warning('Unable to find link url in %s', link)
            continue
        url = match.group(1)
        if url.startswith('/'):
            url = source.asdict()['origin'] + url
        # logging.info('URL %s', url)

        # get the filename for the download file
        file_name_found = False
        save_file_name = None
        attempts = 0
        while attempts < 3:
            response = requests.head(url)
            save_file_name = filename_from_headers(response.headers)
            if save_file_name:
                file_name_found = True
                break
            else:
                logging.info('Failed to find filename via head request, retrying ...')
                attempts += 1
                time.sleep(1)
        if not file_name_found:
            logging.warning('Unable to find filename in headers')
            save_file_name = str(uuid.uuid4()) + '.gz'
        save_file_full = f'{save_file_path}/{save_file_name}'

        # download and save the file
        os.makedirs(save_file_path, exist_ok=True)
        if os.path.isfile(save_file_full):
            logging.info('Skipping existing file: %s', save_file_full)
        else:
            logging.info('Downloading to: %s', save_file_full)
            headers = download_file(url, save_file_full)
            if not file_name_found:
                # try one last time to find filename in headers and rename downloaded file if found
                save_file_name = filename_from_headers(headers)
                if save_file_name:
                    new_file_full = f'{save_file_path}/{save_file_name}'
                    if not os.path.isfile(new_file_full):
                        os.rename(save_file_full, new_file_full)
                        logging.info('Renaming %s to %s', save_file_full, new_file_full)
    return save_file_path


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
