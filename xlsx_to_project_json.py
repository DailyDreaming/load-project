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
    if data.get("contributors.name", []):
        template['contributors'] = []
    for i in range(len(data.get("contributors.name", []))):
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
    if data.get("publications.title", []):
        template['publications'] = []
    for i in range(len(data.get("publications.title", []))):
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
    if data.get("funders.grant_title", []):
        template['funders'] = []
    for i in range(len(data.get("funders.grant_title", []))):
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


def create_cell_suspension_json(data):
    template = {
        "describedBy": "https://schema.humancellatlas.org/type/biomaterial/13.1.1/cell_suspension",
        "schema_type": "biomaterial",
        "estimated_cell_count": round(data['estimated_cell_count'][0]),
        "genus_species": [
            {
                "text": data['selected_cell_type.text'][0],
                "ontology_label": data["genus_species.ontology_label"][0],
                "ontology": data["genus_species.ontology"][0],
            }
        ]
    }
    return template


def process_section(section, prefix='', verify=False):
    project_data = {}
    for i, row in enumerate(section.rows):
        if i == 3:
            cells = [c.value[len(prefix):] if c.value.startswith(prefix) else c.value for c in row if c.value]
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


def parse_project_data_from_xlsx(wb):
    data = process_section(section=wb['Project'], prefix='project.')
    data.update(process_section(section=wb['Project - Publications'], prefix='project.'))
    data.update(process_section(section=wb['Project - Funders'], prefix='project.'))
    data.update(process_section(section=wb['Project - Contributors'], prefix='project.'))
    return data


def parse_cell_suspension_data_from_xlsx(wb):
    data = process_section(section=wb['Cell suspension'], prefix='cell_suspension.')
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
    wb = load_workbook(filename=args.xlsx)
    project_data = parse_project_data_from_xlsx(wb)
    project_json = create_project_json(project_data, uuid=args.uuid, version=timestamp())

    cell_suspension_data = parse_cell_suspension_data_from_xlsx(wb)
    cell_suspension_json = create_cell_suspension_json(cell_suspension_data)

    with open('bundle/project_0.json', 'w') as f:
        f.write(json.dumps(project_json, indent=4))
    print('"bundle/project_0.json" successfully written.')
    with open('bundle/cell_suspension_0.json', 'w') as f:
        f.write(json.dumps(cell_suspension_json, indent=4))
    print('"bundle/cell_suspension_0.json" successfully written.')


if __name__ == "__main__":
    main()
