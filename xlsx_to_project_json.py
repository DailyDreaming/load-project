import json
import sys
import argparse
from datetime import datetime

from openpyxl import load_workbook


def timestamp():
    return datetime.utcnow().strftime("%Y-%m-%dT%H%M%S.%fZ")


def create_project_json(data, uuid, version, verify=False):
    template = {
      "describedBy": "https://schema.humancellatlas.org/type/project/14.1.0/project",
      "schema_type": "project",
      "project_core": {
        "project_short_name": data["project_core.project_short_name"][0],
        "project_title": data["project_core.project_title"][0],
        "project_description": data["project_core.project_description"][0]
      }
    }
    #### Links and Accessions
    optional_includes = ["supplementary_links",
                         "insdc_project_accessions",
                         "geo_series_accessions",
                         "insdc_study_accessions",
                         "biostudies_accessions"]
    for field in optional_includes:
        if data[field][0]:
            template[field] = data[field]

    #### Contributors
    contributors = [i for i in data.get("contributors.name", []) if i]
    if contributors:
        template['contributors'] = []
    for i in range(len(contributors)):
        persons_role = {}
        optional_includes = ["contributors.project_role.text",
                             "contributors.project_role.ontology",
                             "contributors.project_role.ontology_label"]
        for field in optional_includes:
            if data[field][i]:
                persons_role[field.split('.')[-1]] = data[field][i]

        person = {"name": data["contributors.name"][i],
                  "project_role": persons_role}
        optional_includes = ["contributors.email",
                             "contributors.phone",
                             "contributors.institution",
                             "contributors.laboratory",
                             "contributors.address",
                             "contributors.country",
                             "contributors.corresponding_contributor",
                             "contributors.orcid_id"]
        for field in optional_includes:
            if data[field][i]:
                if field == "contributors.corresponding_contributor":
                    if data[field][i] == 'yes':
                        person[field.split('.')[-1]] = True
                    elif data[field][i] == 'no':
                        person[field.split('.')[-1]] = False
                    else:
                        raise RuntimeError('This should never happen.')
                else:
                    person[field.split('.')[-1]] = data[field][i]

        template['contributors'].append(person)

    #### Publications
    publications = [i for i in data.get("publications.title", []) if i]
    if publications:
        template['publications'] = []
    for i in range(len(publications)):
        publication = {"authors": data["publications.authors"][i].split('||'),
                       "title": data["publications.title"][i]}
        optional_includes = ["publications.doi",
                             "publications.pmid",
                             "publications.url"]
        for field in optional_includes:
            if data[field][i]:
                if field == "publications.pmid":
                    publication[field.split('.')[-1]] = int(data[field][i])
                else:
                    publication[field.split('.')[-1]] = data[field][i]
        template['publications'].append(publication)

    #### Funders
    grant_ids = [i for i in data.get("funders.grant_id", []) if i]
    if grant_ids:
        template['funders'] = []
    for i in range(len(grant_ids)):
        funder = {"grant_title": data["funders.grant_title"][i],
                  "grant_id": data["funders.grant_id"][i],
                  "organization": data["funders.organization"][i]}
        template['funders'].append(funder)

    # TODO:  Determine if uuid is made special (from project name combo?).
    # TODO:  Check for previous version.
    template["provenance"] = {
        "document_id": uuid,
        "submission_date": version,
        "update_date": version,
        "schema_major_version": 14,
        "schema_minor_version": 1
        }
    if verify:
        print(json.dumps(template, indent=4))
    return template


def process_section(section, project=True, verify=False):
    project_data = {}
    for i, row in enumerate(section.rows):
        if i == 3:
            cells = [c.value[len('project.'):] for c in row if c.value] if project else [c.value for c in row if c.value]
        if i == 4:
            assert row[0].value == 'FILL OUT INFORMATION BELOW THIS ROW', cells[0]
        if i == 5:
            for j, field in enumerate(cells):
                project_data[field] = [row[j].value]
        if i > 5:
            for j, field in enumerate(cells):
                project_data[field].append(row[j].value)
    if verify:
        print(json.dumps(project_data, indent=4))
    return project_data


def parse_project_data_from_xlsx(file):
    wb = load_workbook(file)
    data = process_section(section=wb['Project'], project=True)
    for project_key in ['Project - Publications', 'Project - Funders', 'Project - Contributors']:
        try:
            data.update(process_section(section=wb[project_key], project=True))
        except KeyError:
            pass
    return data


def main(argv=sys.argv[1:]):
    parser = argparse.ArgumentParser(description='Turn an xlsx file into a project.json file.')
    parser.add_argument("--uuid", type=str,
                        default="4d6f6c96-2a83-43d8-8fe1-0f53bffd4674",  # TODO: Delete this default (Liver Project).
                        help="The project UUID.")
    parser.add_argument("--xlsx", type=str,
                        default='Gary_Bader_9_16.xlsx',  # TODO: Delete this default (Liver Project).
                        help="Path to an xlsx (excel) file.")

    args = parser.parse_args(argv)

    data = parse_project_data_from_xlsx(file='/home/quokka/Downloads/GEOD-93593_HCA_Ontologies_July_2.xlsx')
    project_json = create_project_json(data, uuid=args.uuid, version=timestamp())

    with open('bundle/project_0.json', 'w') as f:
        f.write(json.dumps(project_json, indent=4))
    print('"bundle/project_0.json" successfully written.')


if __name__ == "__main__":
    main()
