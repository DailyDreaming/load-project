## Prod Project Rereival

import requests, pprint, hca, csv, os

azul_project_link = 'https://service.explore.data.humancellatlas.org/repository/projects?size=200'
dss_client = hca.dss.DSSClient(swagger_url='https://dss.data.humancellatlas.org/v1/swagger.json')


def get_bundles(project_uuid: str):
    query = { "query": { "bool": { "must": [ { "match": { "files.project_json.provenance.document_id": project_uuid } } ] } } }
    resp = dss_client.post_search(replica='aws', es_query=query)
    assert resp.get('total_hits') != 0
    bundles = [x.get('bundle_fqid') for x in resp.get('results')]
    return bundles

def get_project_json_uuid(bundle_uuid: str):
    manifest = dss_client.get_bundle(uuid=bundle_uuid,replica='aws')
    files = manifest['bundle']['files']
    found = None
    for file in files:
        if 'project' in file.get('name'):
            found = file
    return found

def get_date_uploaded(file_uuid: str):
    file_metadata = dss_client.get_file(uuid=file_uuid,replica='aws')
    return file_metadata['provenance']['submission_date']

def output_csv(project_list:list):
    file_path = os.path.dirname(__file__)
    with open(f'{file_path}/project_stats.csv', 'w') as f:
        w = csv.DictWriter(f, project_list[1].keys())
        w.writeheader()
        for project in project_list:
            w.writerow(project)

res = requests.get(azul_project_link)
project_info = []

for x in res.json().get('hits'):
    project_uuid = x.get('entryId')
    cellSuspension = [suspension.get('totalCells') for suspension in x.get('cellSuspensions')]
    cell_count = sum(cellSuspension)
    project_info.append({ "project_uuid": project_uuid, "total_cells": cell_count })

    
# need to get ingestion date, can be found from elasticsearch?...  if not we do the long way....
for project in project_info:
    project_bundles = get_bundles(project['project_uuid'])
    file = None
    for fqid in project_bundles:
        bundle,version = fqid.split('.',1)
        file = get_project_json_uuid(bundle)
        if file is not None: #we found a project_json
            break
    date_uploaded = get_date_uploaded(file.get('uuid'))
    project['date_uploaded'] = date_uploaded
    
pprint.pprint(project_info)
output_csv(project_info)
