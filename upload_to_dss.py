import json

from hca import HCAConfig
from hca.dss import DSSClient


hca_config = HCAConfig()

hca_config["DSSClient"].swagger_url = f"https://dss.dev.data.humancellatlas.org/v1/swagger.json"
dss = DSSClient(config=hca_config)

response = dss.upload(src_dir='bundle', replica='aws', staging_bucket='lon-test-data')
print(f'Successful upload.  Bundle information is:\n{json.dumps(response, indent=4)}')
