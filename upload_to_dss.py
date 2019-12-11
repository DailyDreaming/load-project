from uuid import uuid4

from hca import HCAConfig
from hca.dss import DSSClient

hca_config = HCAConfig()

hca_config["DSSClient"].swagger_url = f"https://dss.dev.data.humancellatlas.org/v1/swagger.json"
dss = DSSClient(config=hca_config)

version = dss.create_version()

# TODO: Add Matrix File or... some kind of data.

projects_json_uuid = str(uuid4())
print(f'Uploading project_0.json: {projects_json_uuid}.{version}')
dss.put_file(
    uuid=projects_json_uuid,
    version=version,
    creator_uid=0,
    source_url="s3://lon-test/project_0.json",
)

links_json_uuid = str(uuid4())
print(f'Uploading links.json: {links_json_uuid}.{version}')
dss.put_file(
    uuid=links_json_uuid,
    version=version,
    creator_uid=0,
    source_url="s3://lon-test/links.json",
)

bundle_json_uuid = str(uuid4())
print(f'Uploading bundle containing project_0.json and links.json: {bundle_json_uuid}.{version}')
dss.put_bundle(
    creator_uid=0,
    uuid=bundle_json_uuid,
    version=version,
    replica="aws",
    files=[
        {
            "uuid": projects_json_uuid,
            "version": version,
            "name": "project_0.json",
            "indexed": False,
        },
        {
            "uuid": links_json_uuid,
            "version": version,
            "name": "links.json",
            "indexed": False,
        }
    ],
)

print('Success.')
