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
from copy_static_project import link_project_metadata
from create_project import (
    run,
)
from util import get_target_spreadsheets


xlsxs = get_target_spreadsheets()

for i, xlsx in enumerate(xlsxs):
    print(f'\n% Progress: {i + 1}/{len(xlsxs)} projects ({xlsx.name}).\n'
          f'===========================================================')
    run(xlsx=str(xlsx))

# This will work with:
# for root, dirs, files in os.walk('/home/quokka/yeah/load-project/projects'):
#     for d in dirs:
#         if len(d) == len('4a95101c-9ffc-4f30-a809-f04518a23803'):
# But I think the explicit list is clearer and strictly avoids possible collisions

# TODO: These have no excel sheets... though in theory Will has them and they should be used instead.
static_prod_bundles = [
    '1defdada-a365-44ad-9b29-443b06bd11d6',
    '4a95101c-9ffc-4f30-a809-f04518a23803',
    '4e6f083b-5b9a-4393-9890-2a83da8188f1',
    '005d611a-14d5-4fbf-846e-571a1f874f70',
    '8c3c290d-dfff-4553-8868-54ce45f4ba7f',
    '008e40e8-66ae-43bb-951c-c073a2fa6774',
    '9c20a245-f2c0-43ae-82c9-2232ec6b594f',
    '027c51c6-0719-469f-a7f5-640fe57cbece',
    '74b6d569-3b11-42ef-b6b1-a0454522b4a0',
    '88ec040b-8705-4f77-8f41-f81e57632f7d',
    '091cf39b-01bc-42e5-9437-f419a66c8a45',
    '577c946d-6de5-4b55-a854-cd3fde40bff2',
    '116965f3-f094-4769-9d28-ae675c1b569c',
    '8185730f-4113-40d3-9cc3-929271784c2b',
    'a9c022b4-c771-4468-b769-cabcf9738de3',
    'a29952d9-925e-40f4-8a1c-274f118f1f51',
    'abe1a013-af7a-45ed-8c26-f3793c24a1f4',
    'ae71be1d-ddd8-4feb-9bed-24c3ddb6e1ad',
    'c4077b3c-5c98-4d26-a614-246d12c2e5d7',
    'cc95ff89-2e68-4a08-a234-480eca21ce79',
    'e0009214-c0a0-4a7b-96e2-d6a83e966ce0',
    'f81efc03-9f56-4354-aabb-6ce819c3d414',
    'f83165c5-e2ea-4d15-a5cf-33f3550bffde'
]

for bundle in static_prod_bundles:
    src_dir = Path('projects') / bundle
    link_project_metadata(str(src_dir), file_pattern='*.json')
    print(f'Hard-linked project ("*.json" only) contents: {bundle}/hca to {bundle}/bundle.')
    link_project_metadata(str(src_dir), file_pattern='*.mtx.zip')
