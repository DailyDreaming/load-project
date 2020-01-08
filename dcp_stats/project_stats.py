## Prod Project Rereival

import requests, pprint, hca, csv, os
urls = {"prod" : {"azul":'https://service.explore.data.humancellatlas.org/repository/projects?size=200',
                    "dss": 'https://dss.data.humancellatlas.org/v1/swagger.json'},
        "dev" : {"azul":'https://service.skunk.dev.explore.data.humancellatlas.org/repository/projects?size=200',
                    "dss": 'https://dss.dev.data.humancellatlas.org/v1/swagger.json'}
        }

manual_counts  = {"116965f3-f094-4769-9d28-ae675c1b569c": 15744,
                  "2043c65a-1cf8-4828-a656-9e247d4e64f1": 1733,
                  "4a95101c-9ffc-4f30-a809-f04518a23803": 267360,
                  "f83165c5-e2ea-4d15-a5cf-33f3550bffde": 546183,
                  "abe1a013-af7a-45ed-8c26-f3793c24a1f4": 341079,
                  "a9c022b4-c771-4468-b769-cabcf9738de3": 40993
                  }
#thanks 2 lon <3


def get_bundles(dss_client, project_uuid: str):
    query = { "query": { "bool": { "must": [ { "match": { "files.project_json.provenance.document_id": project_uuid } } ] } } }
    resp = dss_client.post_search(replica='aws', es_query=query)
    assert resp.get('total_hits') != 0
    bundles = [x.get('bundle_fqid') for x in resp.get('results')]
    return bundles

def get_project_json_uuid(dss_client, bundle_uuid: str):
    manifest = dss_client.get_bundle(uuid=bundle_uuid,replica='aws')
    for file in manifest['bundle']['files']:
        if 'project' in file.get('name'):
            return file

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
        if project_uuid in manual_counts:
            cell_count = manual_counts[project_uuid]
        else:
            cellSuspension = [suspension.get('totalCells') for suspension in x.get('cellSuspensions')]
            cell_count = sum(cellSuspension)
        project_info.append({ "project_uuid": project_uuid, "total_cells": cell_count })


    # this is the long way to do it
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
