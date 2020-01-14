"""
Running this will parse all excel files into bundles output into the folder: 'projects/{project_uuid}/bundle'.

This is currently set to loop over two sets of excel files.

One is the set of 6 original project excel files in "data".  These have project_0.json and cell_suspension_0.json 
files included that were downloaded from dss prod to check against.

The second is the set of 71 raw excel files with no comparison in the dss yet (these need to be checked).  These 
also have slightly different formatting (all seem to be missing a "funders" section... and a few are missing the
"publications" section).

Editing 'upload=False' to 'upload=True' with upload to dss dev if you are credentialed to do so.
Otherwise it will just parse the excel files and generate all of the matrix and json files necessary for upload.
"""
from pathlib import Path

from create_project import (
    run,
)
from util import get_target_spreadsheets


xlsxs = get_target_spreadsheets()

for i, xlsx in enumerate(xlsxs):
    print(f'\n% Progress: {i + 1}/{len(xlsxs)} projects ({xlsx.name}).\n'
          f'===========================================================')
    run(xlsx=xlsx.as_posix())
