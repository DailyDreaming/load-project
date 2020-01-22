#!/home/nadove/PycharmProjects/data-portal-summary-stats/.venv/bin/python
from urllib.error import HTTPError

import argparse
from io import (
    BytesIO
)
import json
import logging
from more_itertools import (
    one
)
import os
from pathlib import Path
import sys
from typing import (
    Dict, 
)
from urllib.request import (
    urlopen
)

import boto3
import colorcet as cc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from util import (
    generate_project_uuid,
    get_target_project_dirs,
)

# Local output
output_dir = Path('tSNE')
cache_dir = output_dir / 'cache'


class TSNE:

    base_url = 'https://www.ebi.ac.uk/gxa/sc/'

    @property
    def tsne_url(self) -> str:
        return f'{self.base_url}/experiments/{self.ax_acc}/results/tsne'

    @property
    def clusters_url(self) -> str:
        return f'{self.base_url}/experiment/{self.ax_acc}/download?fileType=cluster&accessKey='

    @property
    def points_url(self) -> str:
        return f'{self.base_url}/json/experiments/{self.ax_acc}/tsneplot/{self.perplexity}/clusters/k/{self.k}'

    def __init__(self, geo_acc: str, perplexity: int):
        self.geo_acc = geo_acc
        self.ax_acc = f'E-GEOD-{self.geo_acc[3:]}'
        self.uuid = generate_project_uuid(self.geo_acc)
        self.perplexity = perplexity
        self.k = self.get_default_k()
        self.clusters = self.get_clusters()

    def get_default_k(self) -> int:
        """
        Obtain the number of clusters SCXA uses by default when coloring by cluster.
        """
        cluster_info = self._cache(
            'cluster_info.tsv',
            lambda f: pd.read_csv(f, sep='\t'),
            lambda cluster_info, f: cluster_info.to_csv(f, sep='\t'),
            lambda: pd.read_csv(BytesIO(urlopen(self.clusters_url).read()), sep='\t')
        )
        return int(one(cluster_info['K'][cluster_info['sel.K']]))

    def get_clusters(self) -> Dict[str, np.ndarray]:
        """
        Get X- and Y- coordiante vectors for tSNE plot points, grouped by cluster.
        :return: clusters names mapped to 2xN arrays where each column is a point and the rows are 
        X- and Y- coordinates.
        """
        points = self._cache(
            'points.json',
            json.load,
            json.dump,
            lambda: json.load(urlopen(self.points_url))
        )
        return {
            cluster["name"]: np.array([
                (point["x"], point["y"])
                for point
                in cluster["data"]
            ]).T
            for cluster
            in points["series"]
        }

    def make_image(self, colormap, save_format: str, dpi: int) -> None:
        """
        Render tSNE using matplotlib.
        """
        fig = plt.figure(figsize=[6, 6])

        plt.title('Clusters', fontweight='bold')
        plt.margins(tight=True)
        
        # Hide axes
        ax = fig.gca()
        ax.spines['bottom'].set_color('0.75')
        for obj in [ax.get_xaxis(), 
                    ax.get_yaxis(),
                    *[ax.spines[side] for side in ('top', 'left', 'right')]]:
            obj.set_visible(False)

        # Make point size shrink with more points
        n_points = np.hstack(list(self.clusters.values())).shape[1]
        marker_size = 225/(n_points**0.65)

        for color, (cluster_id, coords) in zip(colormap, self.clusters.items()):
            ax.scatter(x=coords[0],
                       y=coords[1], 
                       s=marker_size,
                       # Prevent matplotlib warning about ambiguous color format
                       c=np.array(color)[..., np.newaxis],
                       label=cluster_id)
        
        legend = plt.legend(loc='upper center',
                            bbox_to_anchor=(0.5, 0),
                            frameon=False,
                            ncol=4,
                            columnspacing=2.0,
                            labelspacing=0.5)

        # Make sure legened markers are akways the same size regardless of how the
        # plotted points scale with number of points
        for handle in legend.legendHandles:
            handle.set_sizes([25])

        fig.text(0.5, 0.0025, f'tSNE data imported from {self.tsne_url}', ha='center', va='bottom', fontsize='x-small')
        plt.tight_layout()
        plt.savefig(output_dir / f'{self.geo_acc}.{save_format}', dpi=dpi)

    def upload(self, s3_client, bucket):
        image_file = one(str(p) for p in output_dir.iterdir() if p.name.startswith(self.geo_acc))
        image_format = image_file.rsplit('.', 1)[-1]
        key = f'project-assets/project-stats/{self.uuid}/tsne.{image_format}'
        log.info(f'Uploading {self.geo_acc} as {key}')
        s3_client.upload_file(
            Bucket=bucket,
            Key=key,
            Filename=image_file,
            ExtraArgs={
                'ACL': 'public-read',
                'ContentDisposition': 'inline',
                'ContentType': f'image/{image_format}'
            }
        )

    def _cache(self, name, reader, writer, getter):
        path = cache_dir / '_'.join([self.geo_acc, name])
        try:
            with open(path) as f:
                value = reader(f)
        except FileNotFoundError:
            log.info(f'Retreiving {name} (not found in cache)')
            value = getter()
            with open(path, 'w') as f:
                writer(value, f)
        else:
            log.info(f'Loaded {name} from cache')
        return value


if __name__ == "__main__":
    log = logging.getLogger(__name__)
    logging.basicConfig(stream=sys.stdout,
                        level=logging.INFO)

    parser = argparse.ArgumentParser(description='generate tSNE plots from SCXA data and upload to S3', add_help=True)
    # General behavior
    parser.add_argument('-r', '--render-only', action='store_true')
    parser.add_argument('-u', '--upload-only', action='store_true')
    # Filter specific projects
    parser.add_argument('projects', nargs='*')
    # Image properties
    parser.add_argument('--colormap', type=cc.__getattribute__, default='glasbey')
    parser.add_argument('-f', '--image-format', default='png')
    parser.add_argument('--dpi', type=int, default=100)
    # tSNE parameter. Default used by SCXA. Used to generate URLs.
    parser.add_argument('-p', '--perplexity', type=int, default=25)
    # Upload properties
    parser.add_argument('--bucket', default='ux-dev.project-assets.data.humancellatlas.org')
    
    args = parser.parse_args()

    do_render = args.render_only or not args.upload_only
    do_upload = args.upload_only_only or not not args.render_only

    client = boto3.client('s3')

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(cache_dir, exist_ok=True)

    for project_dir in get_target_project_dirs():
        try:
            tsne = TSNE(project_dir.name, args.perplexity)
        except HTTPError:
            log.info(f'Failed to retrieve tSNE data from SCXA for project {project_dir.name}')
        else:
            if do_render:
                tsne.make_image(args.colormap, args.image_format, args.dpi)
            if do_upload:
                tsne.upload(client, args.bucket)
