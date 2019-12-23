import json

import argparse
import sys

from hca import HCAConfig
from hca.dss import DSSClient


def main(argv):
    parser = argparse.ArgumentParser(description='Turn an xlsx file into a project.json file.')
    parser.add_argument("source", type=str,
                        default="bundle",
                        help="The directory to upload as a bundle")

    args = parser.parse_args(argv)
    hca_config = HCAConfig()

    hca_config["DSSClient"].swagger_url = f"https://dss.dev.data.humancellatlas.org/v1/swagger.json"
    dss = DSSClient(config=hca_config)

    response = dss.upload(src_dir=args.source, replica='aws', staging_bucket='lon-test-data')
    print(f'Successful upload.  Bundle information is:\n{json.dumps(response, indent=4)}')


if __name__ == '__main__':
    main(sys.argv[1:])
