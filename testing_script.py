from uuid import uuid4

from xlsx_to_project_json import run

for i in range(6):
    run(uuid=str(uuid4()), xlsx=f'data/test_00{i}.xlsx')
