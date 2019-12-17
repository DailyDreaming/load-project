import uuid
from xlsx_to_project_json import run

# This is used to consistently generate project UUIDs
namespace_uuid = '0887968d-72ec-4c58-bd99-be55953aa462'

for i in range(6):
    run(namespace_uuid, xlsx=f'data/test_00{i}.xlsx', output_dir=f'testing_comparison/test_00{i}')
