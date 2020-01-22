import json
import logging
import re
import unicodedata

import boto3

from _pathlib import Path
from util import get_target_project_dirs

log = logging.getLogger(__name__)


class Main:

    def __init__(self):
        self.s3 = boto3.resource('s3')

    bucket_name = 'ux-dev.project-assets.data.humancellatlas.org'

    key_prefix = 'project-assets/project-matrices/'

    object_acl = {
        'Owner': {
            'DisplayName': 'czi-aws-admins+humancellatlas',
            'ID': '76fe35006be54bbb55cf619bf94684704f14362141f2422a19c3af23e080a148'
        },
        'Grants': [
            {
                'Grantee': {
                    'URI': 'http://acs.amazonaws.com/groups/global/AllUsers',
                    'Type': 'Group'
                },
                'Permission': 'READ'
            }
        ]
    }

    def upload_files_to_bucket(self, project_dir: Path):
        log.info('Checking %s for assets to upload ...', project_dir)
        bundle_dir = project_dir / 'bundle'
        file_extension = '.mtx.zip'
        matrix_file = bundle_dir / ('matrix' + file_extension)
        if matrix_file.exists():
            species_name = self.get_species(bundle_dir)
            if species_name is None:
                log.warning('Could not determine species name, skipping %s.', project_dir)
            else:
                project_uuid = project_dir.name
                file_name = f'{project_uuid}.{species_name}{file_extension}'
                key = self.key_prefix + file_name
                obj = self.s3.Object(self.bucket_name, key)
                try:
                    last_modified = obj.last_modified.timestamp()
                except self.s3.meta.client.exceptions.ClientError as e:
                    if e.response['Error']['Code'] == '404':
                        last_modified = 0
                    else:
                        raise
                # Timezone check unnecessary since S3 date-time is HTTP standards compliant,
                # which uses GMT. Also note Unix time uses UTC as standard.
                if int(matrix_file.stat().st_mtime) >= int(last_modified):
                    content_disposition = f'attachment;filename="{file_name}"'
                    log.info(f'Uploading %s to s3://%s/%s.', matrix_file, self.bucket_name, key)
                    with open(str(matrix_file), 'rb') as f:
                        obj.put(Body=f,
                                ContentDisposition=content_disposition,
                                ContentType='application/zip, application/octet-stream')
                    obj.Acl().put(AccessControlPolicy=self.object_acl)
                else:
                    log.info('%s is up to date.', matrix_file)
        else:
            log.warning('Found no matrix asset for project %s.', project_dir)

    def get_species(self, bundle_dir: Path):
        donor_path = bundle_dir / 'donor_organism_0.json'
        try:
            with open(str(donor_path), 'r') as cs_json:
                cell_suspension_json = json.load(cs_json)
        except FileNotFoundError as e:
            log.warning('Failed to load donor metadata: %s', e)
            return None
        else:
            # FIXME: check other element in genus species array
            try:
                species_name = cell_suspension_json['genus_species'][0]['text']
            except KeyError:
                pass
            except IndexError:
                pass
            else:
                if species_name is not None:
                    species_name = unicodedata.normalize('NFKD', species_name)
                    return re.sub(r'[^\w,.@%&-_()\\[\]/{}]', '_', species_name).strip().lower()
            log.warning('Donor metadata does not specify species name')
            return None

    def run(self):
        for project_dir in get_target_project_dirs(follow_links=True):
            main.upload_files_to_bucket(project_dir)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main = Main()
    main.run()
