from datetime import datetime
import json
import os
import uuid


# TODO: Consolidate similar functions and clean up code.


def timestamp():
    return datetime.utcnow().strftime("%Y-%m-%dT%H%M%S.%fZ")


good_timestamp = timestamp()


def generate_file_uuid(bundle_uuid: str, file_name: str) -> str:
    """
    Deterministically generate a file UUID based on the parent bundle uuid and its file name.
    """
    namespace_uuid = uuid.UUID('4c52e3d0-ffe5-4b4d-a4d0-cb6a6f372b31')
    return str(uuid.uuid5(namespace_uuid, bundle_uuid + file_name))


def get_cell_count(project_uuid):
    with open('/home/quokka/yeah/load-project/.attic/prod_cell_counts.json', 'r') as f:
        counts = json.load(f)
    return counts[project_uuid]['cell_count']


def update_fields(file, filename, bundle_uuid, cell_count, project_uuid=None):
    project_uuid = project_uuid if project_uuid else generate_file_uuid(bundle_uuid, filename)
    with open(file, 'r') as f1:
        update = json.loads(f1.read())

    if filename not in ('links.json', 'cell_suspension_0.json'):
        update['provenance']['document_id'] = project_uuid
        update['provenance']['submission_date'] = good_timestamp
        update['provenance']['update_date'] = good_timestamp

    if filename == 'cell_suspension_0.json':
        update['estimated_cell_count'] = cell_count
        update['describedBy'] = 'https://schema.dev.data.humancellatlas.org/type/biomaterial/13.1.0/cell_suspension'
        update['provenance'] = {'document_id': project_uuid,
                                'submission_date': good_timestamp,
                                'update_date': good_timestamp}

    with open(file, 'w') as f1:
        f1.write(json.dumps(dict(update), indent=4))


def generate_analysis_json(bundle_uuid, output_dir):
    file_name = 'analysis_file_0.json'
    analysis_json = {
        "describedBy": "https://schema.humancellatlas.org/type/file/6.0.0/analysis_file",
        "file_core": {
            "file_name": "matrix.mtx.zip",
            "format": "mtx"
        },
        "schema_type": "file",
        "provenance": {
            "document_id": generate_file_uuid(bundle_uuid, file_name),
            "submission_date": good_timestamp,  # TODO: Fetch from DSS if it exists
            "update_date": good_timestamp
        }
    }
    with open(f'{output_dir}/{file_name}', 'w') as f:
        f.write(json.dumps(analysis_json, indent=4))
    print(f'"{output_dir}/{file_name}" successfully written.')


for root, dirs, _ in os.walk('/home/quokka/yeah/load-project/projects'):
    for d in dirs:
        project_folder = os.path.join(root, d)

        if len(d) == len('4c52e3d0-ffe5-4b4d-a4d0-cb6a6f372b31'):
            cell_count = get_cell_count(d)
            assert cell_count > 1
            for sub_root, _, files in os.walk(project_folder):
                for f in files:
                    file = os.path.join(sub_root, f)
                    if f == 'bundle.json' \
                            or f.startswith('sequence_file') \
                            or (f.startswith('cell_suspension') and f != 'cell_suspension_0.json'):
                        print(f'Deleting: {file}')
                        os.remove(file)
                    elif f == 'links.json':
                        with open(file, 'w') as f1:
                            f1.write(json.dumps({
                                "describedBy": "https://schema.humancellatlas.org/system/1.1.5/links",
                                "schema_type": "link_bundle",
                                "schema_version": "1.1.5",
                                "links": []}, indent=4))
                    elif f.startswith('project_0.json'):
                        print(f'Updating: {file}')
                        update_fields(file=file, filename=f, bundle_uuid=d, cell_count=cell_count, project_uuid=d)
                    elif f.endswith('.json'):
                        print(f'Updating: {file}')
                        print(f'{f}: {cell_count}')
                        update_fields(file=file, filename=f, bundle_uuid=d, cell_count=cell_count)

            matrix_file = os.path.join(project_folder, 'bundle', 'matrix.mtx.zip')
            if os.path.exists(matrix_file):
                generate_analysis_json(bundle_uuid=d, output_dir=os.path.join(project_folder, 'bundle'))
