import json
import logging
from pathlib import Path
from uuid import UUID

from hca import HCAConfig
from hca.dss import DSSClient

log = logging.getLogger(__name__)


def main():
    logging.basicConfig(level=logging.INFO)
    hca_config = HCAConfig()
    hca_config["DSSClient"].swagger_url = f"https://dss.dev.data.humancellatlas.org/v1/swagger.json"
    dss = DSSClient(config=hca_config)

    projects = Path.cwd() / 'projects'
    for project in projects.iterdir():
        if project.is_dir() and not project.is_symlink():
            log.info('Uploading %s', project)
            bundle_uuid = project.name
            assert str(UUID(bundle_uuid)) == bundle_uuid
            bundle = project / 'bundle'
            if bundle.is_dir():
                response = dss.upload(src_dir=bundle.as_posix(),
                                      replica='aws',
                                      staging_bucket='lon-test-data',
                                      bundle_uuid=bundle_uuid)
                print(f'Successful upload.  Bundle information is:\n{json.dumps(response, indent=4)}')
            else:
                log.warning('Skipping %s because metadata is missing', project)


if __name__ == '__main__':
    main()
