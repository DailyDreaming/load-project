## Prod Project Rereival

import requests, pprint, hca, csv, os
urls = { "prod" : {"azul":'https://service.explore.data.humancellatlas.org/repository/projects?size=200',
                    "dss": 'https://dss.data.humancellatlas.org/v1/swagger.json'},
         "dev" : {"azul":'https://service.skunk.dev.explore.data.humancellatlas.org/repository/projects?size=200',
                    "dss": 'https://dss.dev.data.humancellatlas.org/v1/swagger.json'},
         }




def get_bundles(dss_client, project_uuid: str):
    query = { "query": { "bool": { "must": [ { "match": { "files.project_json.provenance.document_id": project_uuid } } ] } } }
    resp = dss_client.post_search(replica='aws', es_query=query)
    assert resp.get('total_hits') != 0
    bundles = [x.get('bundle_fqid') for x in resp.get('results')]
    return bundles

def get_project_json_uuid(dss_client, bundle_uuid: str):
    manifest = dss_client.get_bundle(uuid=bundle_uuid,replica='aws')
    files = manifest['bundle']['files']
    found = None
    for file in files:
        if 'project' in file.get('name'):
            found = file
    return found

def get_date_submission(dss_client, file_uuid: str):
    file_metadata = dss_client.get_file(uuid=file_uuid,replica='aws')
    return file_metadata['provenance']['submission_date']

def output_csv(project_list:list,stage):
    file_path = os.path.dirname(__file__)
    with open(f'{file_path}/project-stats-{stage}.csv', 'w') as f:
        w = csv.DictWriter(f, project_list[1].keys())
        w.writeheader()
        for project in project_list:
            w.writerow(project)


def run(stage):
    res = requests.get(urls[stage]['azul'])
    dss_client = hca.dss.DSSClient(swagger_url=urls[stage]['dss'])
    project_info = []

    # get all the projects
    for x in res.json().get('hits'):
        project_uuid = x.get('entryId')
        cellSuspension = [suspension.get('totalCells') for suspension in x.get('cellSuspensions')]
        cell_count = sum(cellSuspension)
        project_info.append({ "project_uuid": project_uuid, "total_cells": cell_count })

    # need to get ingestion date, can be found from elasticsearch?...  if not we do the long way....
    for project in project_info:
        project_bundles = get_bundles(dss_client, project['project_uuid'])
        file = None
        for fqid in project_bundles:
            bundle,version = fqid.split('.',1)
            file = get_project_json_uuid(dss_client, bundle)
            if file is not None: #we found a project_json
                break
        date_submission = get_date_submission(dss_client, file.get('uuid'))
        project['date_submission'] = date_submission.split('T')[0]

    sorted_project_info = sorted(project_info,key=lambda d: d['date_submission'])
    for idx,project  in enumerate(sorted_project_info):
        project['total_projects'] = idx + 1

    pprint.pprint(sorted_project_info)
    output_csv(sorted_project_info,stage)

if __name__ == "__main__":
    run('prod')
    run('dev')