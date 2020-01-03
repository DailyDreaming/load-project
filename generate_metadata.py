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
    get_spreadsheet_paths,
    run,
)

for sub_dir in 'existing', 'new':
    src_dir = Path('spreadsheets') / sub_dir
    projects = get_spreadsheet_paths(src_dir)
    for i, project in enumerate(projects):
        print(f'\n% Progress: {i + 1}/{len(projects)} projects ({project}).\n'
              f'===========================================================')
        run(xlsx=project.as_posix())
