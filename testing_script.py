from xlsx_to_project_json import run, namespace_uuid

for i in range(6):
    run(namespace_uuid, xlsx=f'data/test_00{i}.xlsx', output_dir=f'testing_comparison/test_00{i}')
