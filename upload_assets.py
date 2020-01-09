import argparse
import json
import logging
import os
import re
import sys
import unicodedata

import boto3

log = logging.getLogger(__name__)


class ProjectMatrixUploader:

    def __init__(self):
        self.s3 = boto3.resource('s3')

    def upload_files_to_bucket(self, bucket_name, geo_projects_directory):
        key_prefix = 'project-assets/project-matrices/'
        for project_uuid in os.listdir(geo_projects_directory):
            project_dir = os.path.join(geo_projects_directory, project_uuid)
            bundle_dir = os.path.join(project_dir, 'bundle')
            # FIXME add later to check other donor_organisim files
            file_extension = '.mtx.zip'
            matrix_path = os.path.join(bundle_dir, 'matrix' + file_extension)
            if not os.path.islink(project_dir):
                if os.path.exists(matrix_path):
                    species_name = self.get_species(bundle_dir)
                    if species_name is None:
                        log.warning('Could not determine species name. Skipping project %s.', project_uuid)
                    else:
                        key = f'{key_prefix}{project_uuid}.{species_name}{file_extension}'
                        content_disposition = f'attachment;filename="{project_uuid}.{species_name}{file_extension}"'
                        obj = self.s3.Object(bucket_name, key)

                        log.info(f'Uploading %s to s3://%s/%s.', matrix_path, bucket_name, key)
                        with open(matrix_path, 'rb') as mf:
                            obj.put(Body=mf,
                                    ContentDisposition=content_disposition,
                                    ContentType='application/zip, application/octet-stream')
                        obj.Acl().put(AccessControlPolicy={
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
                        })
                else:
                    log.warning('Found no matrix asset for project %s.', project_uuid)

    def get_species(self, bundle_dir):
        donor_path = os.path.join(bundle_dir, 'donor_organism_0.json')
        try:
            with open(donor_path, 'r') as cs_json:
                cell_suspension_json = json.load(cs_json)
                # FIXME check other element in genus species array
                species_name = cell_suspension_json['genus_species'][0]['text']
                species_name = unicodedata.normalize('NFKD', species_name)
                return re.sub(r'[^\w,.@%&-_()\\[\]/{}]', '_', species_name).strip().lower()
        except FileNotFoundError:
            log.warning('Failed to load donor metadata', exc_info=True)
            return None


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser('Uploads all matrix files from a single directory to a S3 bucket.')
    parser.add_argument('-b', '--bucket', default='ux-dev.project-assets.data.humancellatlas.org')
    parser.add_argument('-d', '--directory', default='projects')
    options = parser.parse_args(sys.argv[1:])
    uploader = ProjectMatrixUploader()
    uploader.upload_files_to_bucket(options.bucket, options.directory)
