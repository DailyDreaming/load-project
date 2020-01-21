import argparse
import logging
import shutil
import sys
from typing import (
    Sequence,
)

from _pathlib import Path
from util import get_target_project_dirs

log = logging.getLogger(__name__)


class Main:

    def __init__(self, args: Sequence[str]) -> None:
        super().__init__()
        parser = argparse.ArgumentParser()
        parser.add_argument('--dry-run', '--dryrun', '-n', default=False, action='store_true')
        parser.add_argument('artifacts', metavar='GLOB', nargs="+")
        self.args = parser.parse_args(args)

    def clean_project(self, project_dir: Path):
        log.info('Looking for artifacts to clean in project %s. ...', project_dir)
        for glob in self.args.artifacts:
            for artifact in project_dir.glob(glob):
                if artifact.is_dir():
                    if self.args.dry_run:
                        log.info('    Would recursively remove directory %s', artifact)
                    else:
                        log.info('    Recursively removing directory %s', artifact)
                        shutil.rmtree(artifact)
                else:
                    if self.args.dry_run:
                        log.info('    Would remove file %s', artifact)
                    else:
                        log.info('    Removing file %s', artifact)
                        artifact.unlink()

    def run(self):
        projects = get_target_project_dirs()
        for project in projects:
            self.clean_project(project)


if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(levelname)s %(threadName)s: %(message)s',
                        level=logging.INFO)
    Main(sys.argv[1:]).run()
