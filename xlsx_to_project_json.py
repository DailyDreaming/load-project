import json
import os
import sys
import argparse
import uuid
from datetime import datetime

from download_geo_matrix import download_supplementary_files, extract_rename_and_zip

from openpyxl import load_workbook


def timestamp():
    return datetime.utcnow().strftime("%Y-%m-%dT%H%M%S.%fZ")


def create_project_json(data, namespace_uuid, version, verify=False):
    project_json = {
      "describedBy": "https://schema.humancellatlas.org/type/project/14.1.0/project",
      "schema_type": "project",
      "project_core": {
        "project_short_name": data["project_core.project_short_name"][0],
        "project_title": data["project_core.project_title"][0],
        "project_description": data["project_core.project_description"][0]
      }
    }
    # Links and Accessions
    optional_includes = ["supplementary_links",
                         "insdc_project_accessions",
                         "geo_series_accessions",
                         "insdc_study_accessions",
                         "biostudies_accessions",
                         "array_express_accessions"]
    for field in optional_includes:
        if field in data:
            if data[field][0]:
                project_json[field] = data[field][0].split('||')
            if len(data[field]) > 1:
                raise RuntimeError('This should never happen.')

    # Contributors
    contributors = [i for i in data.get("contributors.name", []) if i]
    if contributors:
        project_json['contributors'] = []
    for i in range(len(contributors)):
        persons_role = {}
        optional_includes = ["contributors.project_role.text",
                             "contributors.project_role.ontology",
                             "contributors.project_role.ontology_label"]
        for field in optional_includes:
            if data[field][i]:
                persons_role[field.split('.')[-1]] = data[field][i]

        person = {"name": data["contributors.name"][i]}
        if persons_role:
            person["project_role"] = persons_role
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

        project_json['contributors'].append(person)

    # Publications
    publications = [i for i in data.get("publications.title", []) if i]
    if publications:
        project_json['publications'] = []
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
        project_json['publications'].append(publication)

    # Funders
    grant_ids = [i for i in data.get("funders.grant_id", []) if i]
    if grant_ids:
        project_json['funders'] = []
    funders = []
    for i in range(len(grant_ids)):
        funder = {}
        if data["funders.grant_title"][i]:
            funder['grant_title'] = data["funders.grant_title"][i]
        if data["funders.grant_id"][i]:
            funder['grant_id'] = data["funders.grant_id"][i]
        if data["funders.organization"][i]:
            funder['organization'] = data["funders.organization"][i]
        if funder:
            funders.append(funder)
    if funders:
        project_json['funders'] = funders

    deterministic_uuid = uuid.uuid5(uuid.UUID(namespace_uuid), ''.join(project_json['geo_series_accessions']))
    project_json["provenance"] = {
        "document_id": str(deterministic_uuid),
        "submission_date": version,
        "update_date": version,
        "schema_major_version": 14,
        "schema_minor_version": 1
        }
    if verify:
        print(json.dumps(project_json, indent=4))
    return project_json, deterministic_uuid


def create_cell_suspension_jsons(data):
    version = timestamp()
    cell_suspension_json = {
        "describedBy": "https://schema.humancellatlas.org/type/biomaterial/13.1.1/cell_suspension",
        "schema_type": "biomaterial",
        "estimated_cell_count": round(data.get('estimated_cell_count', [1])[0]),
        "biomaterial_core": {
            "biomaterial_id": data['biomaterial_core.biomaterial_id'][0],
            "biomaterial_description": data['biomaterial_core.biomaterial_description'][0],
            "ncbi_taxon_id": [
                data['biomaterial_core.ncbi_taxon_id'][0]
            ]
        },
        "genus_species": [
            {
                "text": data['genus_species.text'][0],
                "ontology_label": data["genus_species.ontology_label"][0],
                "ontology": data["genus_species.ontology"][0],
            }
        ],
        "provenance": {
            "document_id": str(uuid.uuid4()),
            "submission_date": version,
            "update_date": version,
            "schema_major_version": 13,
            "schema_minor_version": 3
        }
    }
    return cell_suspension_json


def is_known_divider(row_value, right_after_key_declaration):
    if 'BELOW THIS ROW' in row_value or 'below this line' in row_value:
        return True
    if right_after_key_declaration and not row_value:
        return True


def process_section(section, prefix='', verify=False):
    project_data = {}
    header_section = True
    right_after_key_declaration = False
    cells = []
    for i, row in enumerate(section.rows):

        if not header_section:
            for j, field in enumerate(cells):
                if field in project_data:
                    project_data[field].append(row[j].value)
                else:
                    project_data[field] = [row[j].value]

        row_value = row[0].value or ''
        if header_section and row_value.startswith(prefix):
            cells = [c.value[len(prefix):] if c.value.startswith(prefix) else c.value for c in row if c.value]
            right_after_key_declaration = True

        if is_known_divider(row_value, right_after_key_declaration):
            header_section = False
    if verify:
        print(json.dumps(project_data, indent=4))
    return project_data


def parse_project_data_from_xlsx(wb):
    data = process_section(section=wb['Project'], prefix='project.')
    for project_key in ['Project - Publications',
                        'Project - Publication',
                        'Project - Funders',
                        'Project - Funder',
                        'Project - Contributors',
                        'Project - Contributor']:
        try:
            data.update(process_section(section=wb[project_key], prefix='project.'))
            print(f'{project_key} found in the data.')
        except KeyError:
            print(f'{project_key} not found in the data!')
    return dict(data)


def parse_cell_suspension_data_from_xlsx(wb):
    data = process_section(section=wb['Cell suspension'], prefix='cell_suspension.')
    return data


def add_matrix_file(accessions, project_uuid, out_dir):
    for acc in accessions:
        download_dir = download_supplementary_files(acc)
        if download_dir:
            matching_files = [f for f in sorted(os.listdir(download_dir)) if f.startswith(acc) and f.endswith('.gz')]
            if matching_files:
                print(f'Matching files found for {acc}: {matching_files}')
                zip_files = extract_rename_and_zip(download_dir, matching_files, project_uuid)
                if zip_files:
                    # move the zip file to the correct location
                    file = zip_files[0]
                    if not os.path.isfile(f'{out_dir}/{file}'):
                        os.rename(f'{download_dir}/{file}', f'{out_dir}/{file}')
                break
        print(f'No matching files found for {acc}')


def generate_project_json(wb, namepace_uuid, output_dir):
    project_data = parse_project_data_from_xlsx(wb)
    project_json, project_uuid = create_project_json(project_data, namespace_uuid=namepace_uuid, version=timestamp())
    with open(f'{output_dir}/project_0.json', 'w') as f:
        f.write(json.dumps(project_json, indent=4))
    print(f'"{output_dir}/project_0.json" successfully written.')
    return project_json, project_uuid


def generate_cell_suspension_json(wb, output_dir):
    cell_suspension_data = parse_cell_suspension_data_from_xlsx(wb)
    cell_json = create_cell_suspension_jsons(cell_suspension_data)
    with open(f'{output_dir}/cell_suspension_0.json', 'w') as f:
        f.write(json.dumps(cell_json, indent=4))
    print(f'"{output_dir}/cell_suspension_0.json" successfully written.')


def generate_links_json(output_dir):
    links = {
        'describedBy': 'https://schema.humancellatlas.org/system/1.1.5/links',
        'schema_type': 'link_bundle',
        'schema_version': '1.1.5',
        'links': []
    }
    with open(f'{output_dir}/links.json', 'w') as f:
        f.write(json.dumps(links, indent=4))
    print(f'"{output_dir}/links.json" successfully written.')


def run(namepace_uuid, xlsx, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    wb = load_workbook(xlsx)

    project_json, project_uuid = generate_project_json(wb, namepace_uuid, output_dir)
    generate_cell_suspension_json(wb, output_dir)
    generate_links_json(output_dir)

    add_matrix_file(project_json['geo_series_accessions'], project_uuid, output_dir)


def main(argv):
    parser = argparse.ArgumentParser(description='Turn an xlsx file into a project.json file.')
    parser.add_argument("--uuid", type=str,
                        default="4d6f6c96-2a83-43d8-8fe1-0f53bffd4674",  # TODO: Delete this.
                        help="The project UUID.")
    parser.add_argument("--xlsx", type=str,
                        default='data/test_project_000.xlsx',  # TODO: Delete this.
                        help="Path to an xlsx (excel) file.")
    parser.add_argument("--output_dir", type=str,
                        default='bundle',
                        help="Path to an output directory.")

    args = parser.parse_args(argv)
    run(args.uuid, args.xlsx, args.output_dir)


if __name__ == "__main__":
    main(sys.argv[1:])
