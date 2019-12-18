import os

from uuid import uuid4
from xlsx_to_project_json import run


# # This script only tests parsing (NO uploading done) so the uuid doesn't matter.
# namespace_uuid = '0887968d-72ec-4c58-bd99-be55953aa462'
#
# for i in range(6):
#     print(f'\nProject: {i}')
#     run(str(uuid4()), xlsx=f'data/test_00{i}.xlsx', output_dir=f'testing_comparison/test_00{i}')

mypath = 'raw_excel_inputs'
j = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f)) and f.endswith('.0.xlsx')]

for i in j:
    print(f'\nProject: {i}')
    run(str(uuid4()), xlsx=f'raw_excel_inputs/{i}', output_dir=f'testing_comparison/{i[:-5]}', upload=True)