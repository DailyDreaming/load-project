"""
Running this will parse all excel files into bundles output into the folder: 'testing_comparison'.

This is currently set to loop over two sets of excel files.

One is the set of 6 original project excel files in "data".  These have project_0.json and cell_suspension_0.json 
files included that were downloaded from dss prod to check against.

The second is the set of 71 raw excel files with no comparison in the dss yet (these need to be checked).  These 
also have slightly different formatting (all seem to be missing a "funders" section... and a few are missing the
"publications" section).

Editing 'upload=False' to 'upload=True' with upload to dss dev if you are credentialed to do so.
Otherwise it will just parse the excel files and generate all of the matrix and json files necessary for upload.
"""
import os
import shutil

from create_project import run

if os.path.exists('testing_comparison/'):
    shutil.rmtree('testing_comparison/')

"""6 ORIGINAL DATASETS (ALREADY IN THE DSS)"""
# These are the original excel files provided that currently exist in dss prod
# and we have finished examples to compare against.
for project in range(6):
    print(f'\nProject: test_00{project}')
    run(xlsx=f'data/test_00{project}.xlsx',
        # output_dir=f'testing_comparison/test_00{project}',
        upload=False)

"""71 RAW DATASETS (STATUS NOT PARSED)"""
# Downloaded from a spreadsheet of spreadsheets and assumed to be (mostly) complete projects.
# These inputs were provided with the labels "finished" or "full".
# Differences assumed are inferred from skimming over the files.
# I chose to use the inputs which end in ".0.xlsx" ("finished") rather than the normal ".xlsx" extension ("full").
#
# These are missing fields such as the "funders" section (as opposed to the 6 excel files above).
# Not sure of other differences yet.
src_dir = 'raw_excel_inputs'
projects = [f for f in os.listdir(src_dir) if os.path.isfile(os.path.join(src_dir, f)) and f.endswith('.0.xlsx')]

for i, project in enumerate(projects):
    print(f'\n% Progress: {i + 1}/{len(projects)} projects ({project}).\n'
          f'===========================================================')
    run(xlsx=f'raw_excel_inputs/{project}',
        # output_dir=f'testing_comparison/{project[:-5]}',
        upload=False)
