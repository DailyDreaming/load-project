import argparse
import boto3
import json
import re
import os
import logging
import unicodedata
import sys

log = logging.getLogger(__name__)


class ProjectMatrixUploader:
    def __init__(self):
        self.s3 = boto3.resource('s3')

    def upload_files_to_bucket(self, bucket_name, geo_projects_directory):
        key_prefix = 'project-assets/project-matrices/'
        for project_uuid in os.listdir(geo_projects_directory):
            project_directory = os.path.join(geo_projects_directory, project_uuid)
            bundle_directory = os.path.join(project_directory, 'bundle')
            # FIXME add later to check other donor_organisim files
            cs_json_path = os.path.join(bundle_directory, 'donor_organism_0.json')
            file_extension = '.mtx.zip'
            matrix_path = os.path.join(bundle_directory, 'matrix' + file_extension)
            if not os.path.islink(project_directory) and os.path.exists(matrix_path):
                with open(cs_json_path, 'r') as cs_json:
                    cell_suspension_json = json.load(cs_json)
                    # FIXME check other element in genus species array
                    species_name = cell_suspension_json['genus_species'][0]['text']
                    species_name = unicodedata.normalize('NFKD', species_name)
                    species_name = re.sub(r'[^\w,.@%&-_()\\[\]/{}]', '_', species_name).strip().lower()
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


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser('Uploads all matrix files from a single directory to a S3 bucket.')
    parser.add_argument('-b', '--bucket', default='ux-dev.project-assets.data.humancellatlas.org')
    parser.add_argument('-d', '--directory', default='projects')
    options = parser.parse_args(sys.argv[1:])
    uploader = ProjectMatrixUploader()
    uploader.upload_files_to_bucket(options.bucket, options.directory)
