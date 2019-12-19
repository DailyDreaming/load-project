import json
import os
import sys
import argparse
import copy
from typing import Sequence
import uuid
from datetime import datetime

from hca import HCAConfig
from hca.dss import DSSClient

from download_geo_matrix import (
    download_supplementary_files,
    extract_rename_and_zip,
    extract_tar_files
)

from geo_namespace import deterministic_uuid

from openpyxl import load_workbook


# TODO: Consolidate similar functions and clean up code.
# TODO: Make file uuid deterministic based on bundle uuid and file name


def timestamp():
    return datetime.utcnow().strftime("%Y-%m-%dT%H%M%S.%fZ")


def generate_project_uuid(geo_accessions: Sequence[str]) -> str:
    """
    Deterministically generate a project UUID based on a series of GEO accession
    numbers.
    """
    # Chosen to match https://github.com/DailyDreaming/load-project/pull/6/files#diff-909893068c9afb3ccbb8ea268955eb7dR30
    # at the suggestion of Daniel Sotirhos.
    # It's essential that this is a hardcoded constant so that other scripts can
    # deterministically derive the same HCA UUIDs from GEO accessions.
    return deterministic_uuid(''.join(sorted(geo_accessions)))


def generate_file_uuid(bundle_uuid: str, file_name: str) -> str:
    """Deterministically generate a file UUID based on the parent bundle uuid and its file name."""
    namespace_uuid = uuid.UUID(bundle_uuid)
    return str(uuid.uuid5(namespace_uuid, file_name))


def get_cell_counts():
    with open('cell_counts.json', 'r') as f:
        cell_counts = json.loads(f.read())
    return cell_counts


def create_project_json(data, version, verify=False):
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

    project_uuid = generate_project_uuid(project_json['geo_series_accessions'])
    project_json["provenance"] = {
        "document_id": str(project_uuid),
        "submission_date": version,
        "update_date": version,
        "schema_major_version": 14,
        "schema_minor_version": 1
        }
    if verify:
        print(json.dumps(project_json, indent=4))
    return project_json, project_uuid


def create_cell_suspension_jsons(data, cell_count, file_uuid, i=0):
    version = timestamp()
    cell_suspension_json = {
        "describedBy": "https://schema.humancellatlas.org/type/biomaterial/13.1.0/cell_suspension",
        "schema_type": "biomaterial"
    }

    cell_suspension_json.update(fill_sections(data=data,
                                              as_list=['genus_species',
                                                       'selected_cell_types'],
                                              keys={
                                                  'biomaterial_core':
                                                      ['biomaterial_id',
                                                       'biomaterial_name',
                                                       'biomaterial_description',
                                                       'ncbi_taxon_id',
                                                       'biosamples_accession',
                                                       'insdc_sample_accession',
                                                       'genotype'],
                                                  'genus_species':
                                                      ['text',
                                                       'ontology',
                                                       'ontology_label'],
                                                  'cell_morphology':
                                                      ['percent_cell_viability',
                                                       'cell_viability_method'],
                                                  'selected_cell_types':
                                                      ['text',
                                                       'ontology',
                                                       'ontology_label'],
                                              },
                                              i=i))
    cell_suspension_json['estimated_cell_count'] = cell_count

    cell_suspension_json["provenance"] = {
        "document_id": str(file_uuid),
        "submission_date": version,  # TODO: Fetch from DSS if it exists
        "update_date": version
    }
    return cell_suspension_json


def fill_sections(data, keys, i, as_list=[]):
    # TODO: Add types
    # TODO: Handle lists intelligently
    document = {}
    for primary_key, secondary_keys in keys.items():
        if secondary_keys:
            document.update(conditionally_create_subdict(data,
                                                         primary_key=primary_key,
                                                         secondary_keys=secondary_keys,
                                                         i=i,
                                                         use_list=bool(primary_key in as_list)))
        else:
            if primary_key in data:
                if data[primary_key][i]:
                    if primary_key == 'organism_age':
                        document[primary_key] = str(data[primary_key][i])
                    else:
                        document[primary_key] = data[primary_key][i]
    return document


def conditionally_create_subdict(data, primary_key, secondary_keys, i, use_list=False):
    prime_dictionary = {}
    sub_dictionary = {}
    for key in secondary_keys:
        if f'{primary_key}.{key}' in data:
            if len([j for j in data[f'{primary_key}.{key}'] if j]) > 0:
                if key == 'ncbi_taxon_id':
                    sub_dictionary[key] = parse_ncbi_taxon_ids(data[f'{primary_key}.{key}'][i])
                else:
                    sub_dictionary[key] = data[f'{primary_key}.{key}'][i]
    if sub_dictionary:
        prime_dictionary[primary_key.split('.')[-1]] = [sub_dictionary] if use_list else sub_dictionary
    return prime_dictionary


def create_specimen_from_organism_json(data, file_uuid, i=0):
    version = timestamp()
    specimen_from_organism_json = {
        "describedBy": "https://schema.humancellatlas.org/type/biomaterial/10.2.0/specimen_from_organism",
        "schema_type": "biomaterial"
    }
    specimen_from_organism_json.update(
        fill_sections(
            data=data,
            as_list=[
                'genus_species',
                'organ_parts',
                'diseases'
            ],
            keys={
                'biomaterial_core':
                    [
                        'biomaterial_id',
                        'biomaterial_name',
                        'biomaterial_description',
                        'ncbi_taxon_id',
                        'biosamples_accession',
                        'insdc_sample_accession',
                        'genotype'
                    ],
                'genus_species':
                    [
                        'text',
                        'ontology',
                        'ontology_label'
                    ],
                'organ':
                    [
                        'text',
                        'ontology',
                        'ontology_label'
                    ],
                'organ_parts':
                    [
                        'text',
                        'ontology',
                        'ontology_label'
                    ],
                'diseases':
                    [
                        'text',
                        'ontology',
                        'ontology_label'
                    ],
                'preservation_storage':
                    [
                        'preservation_method',
                        'storage_method'
                    ]
            },
            i=i
        )
    )
    specimen_from_organism_json["provenance"] = {
            "document_id": file_uuid,
            "submission_date": version,  # TODO: Fetch from DSS if it exists
            "update_date": version
        }
    return specimen_from_organism_json


def parse_ncbi_taxon_ids(ids):
    if isinstance(ids, int):
        return [ids]
    elif isinstance(ids, float):
        return [int(ids)]
    elif isinstance(ids, str):
        return [int(i) for i in ids.split(',') if i.strip()]
    else:
        raise RuntimeError('This should never happen.')


def create_donor_organism_json(data, file_uuid, i=0):
    version = timestamp()
    donor_organism_json = {
        "describedBy": "https://schema.humancellatlas.org/type/biomaterial/15.3.0/donor_organism",
        "schema_type": "biomaterial"
    }
    donor_organism_json.update(
        fill_sections(
            data=data,
            as_list=[
                'genus_species',
                'diseases'
            ],
            keys={
                'biomaterial_core':
                    [
                        'biomaterial_id',
                        'biomaterial_name',
                        'biomaterial_description',
                        'ncbi_taxon_id',
                        'biosamples_accession',
                        'insdc_sample_accession',
                        'genotype'],
                'genus_species':
                    [
                        'text',
                        'ontology',
                        'ontology_label'
                    ],
                'is_living': [],
                'sex': [],
                'diseases':
                    [
                        'text',
                        'ontology',
                        'ontology_label'
                    ],
                'organism_age': [],
                'organism_age_unit':
                    [
                        'text',
                        'ontology',
                        'ontology_label'
                    ],
                'human_specific':
                    [
                        'body_mass_index'
                    ],
                'death':
                    [
                        'cause_of_death'
                    ],
                'development_stage':
                    [
                        'text',
                        'ontology',
                        'ontology_label'
                    ]
            },
            i=i
        )
    )
    ethnicity = conditionally_create_subdict(
        data=data,
        primary_key='human_specific.ethnicity',
        secondary_keys=[
            'text',
            'ontology',
            'ontology_label'
        ],
        i=i,
        use_list=True
    )
    if ethnicity:
        donor_organism_json.get('human_specific', {}).update(ethnicity)
    donor_organism_json["provenance"] = {
            "document_id": file_uuid,
            "submission_date": version,  # TODO: Fetch from DSS if it exists
            "update_date": version
        }
    return donor_organism_json


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


def get_harmonized_project_sections(wb, section_keywords):
    """The excel fields change randomly and very slightly for no reason so we only look for key words."""
    project_data = {}
    for key in wb:
        for section in section_keywords:
            if section in key.title.lower():
                project_data[section] = wb[str(key.title)]
    return project_data


def parse_project_data_from_xlsx(wb):
    data = process_section(section=wb['Project'], prefix='project.')

    section_keywords = ['publication', 'funder', 'contributor']
    project_sections = get_harmonized_project_sections(wb, section_keywords)
    for project_key in section_keywords:
        try:
            data.update(process_section(section=project_sections[project_key], prefix='project.'))
        except KeyError:
            print(f'{project_key} section not found in the data!')
    return data


def parse_cell_suspension_data_from_xlsx(wb):
    data = {}

    section_keywords = ['suspension']
    project_sections = get_harmonized_project_sections(wb, section_keywords)
    for project_key in section_keywords:
        try:
            data.update(process_section(section=project_sections[project_key], prefix='cell_suspension.'))
        except KeyError:
            print(f'{project_key} section not found in the data!')
    return data


def parse_specimen_from_organism_data_from_xlsx(wb):
    data = {}

    section_keywords = ['specimen']
    project_sections = get_harmonized_project_sections(wb, section_keywords)
    for project_key in section_keywords:
        try:
            data.update(process_section(section=project_sections[project_key], prefix='specimen_from_organism.'))
        except KeyError:
            print(f'{project_key} section not found in the data!')
    return data


def parse_donor_organism_data_from_xlsx(wb):
    data = {}

    section_keywords = ['donor']
    project_sections = get_harmonized_project_sections(wb, section_keywords)
    for project_key in section_keywords:
        try:
            data.update(process_section(section=project_sections[project_key], prefix='donor_organism.'))
        except KeyError:
            print(f'{project_key} section not found in the data!')
    return data


def add_matrix_file(accessions, project_uuid, out_dir):
    for acc in accessions:
        download_dir = download_supplementary_files(acc)
        if download_dir:
            extract_tar_files(download_dir)
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


def generate_project_json(wb, output_dir):
    project_data = parse_project_data_from_xlsx(wb)
    project_json, project_uuid = create_project_json(project_data, version=timestamp())
    with open(f'{output_dir}/project_0.json', 'w') as f:
        f.write(json.dumps(project_json, indent=4))
    print(f'"{output_dir}/project_0.json" successfully written.')
    return project_json, project_uuid


def generate_cell_suspension_json(wb, output_dir, cell_count, bundle_uuid):
    file_name = 'cell_suspension_0.json'
    cell_suspension_data = parse_cell_suspension_data_from_xlsx(wb)
    cell_json = create_cell_suspension_jsons(data=cell_suspension_data,
                                             cell_count=cell_count,
                                             file_uuid=generate_file_uuid(bundle_uuid, file_name))
    with open(f'{output_dir}/{file_name}', 'w') as f:
        f.write(json.dumps(cell_json, indent=4))
    print(f'"{output_dir}/{file_name}" successfully written.')


def generate_specimen_from_organism_json(wb, output_dir, bundle_uuid):
    file_name = 'specimen_from_organism_0.json'
    specimen_from_organism_data = parse_specimen_from_organism_data_from_xlsx(wb)
    specimen_from_organism_json = create_specimen_from_organism_json(
        data=specimen_from_organism_data,
        file_uuid=generate_file_uuid(bundle_uuid, file_name))
    with open(f'{output_dir}/{file_name}', 'w') as f:
        f.write(json.dumps(specimen_from_organism_json, indent=4))
    print(f'"{output_dir}/{file_name}" successfully written.')


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


def generate_donor_organism_jsons(wb, output_dir, bundle_uuid):
    donor_organism_data = parse_donor_organism_data_from_xlsx(wb)
    donors = [donor for donor in donor_organism_data['biomaterial_core.biomaterial_id'] if donor]
    for donor_number in range(len(donors)):
        generate_donor_organism_json(donor_organism_data, output_dir, donor_number, bundle_uuid)


def generate_donor_organism_json(data, output_dir, donor_number, bundle_uuid):
    file_name = f'donor_organism_{donor_number}.json'
    donor_organism_json = create_donor_organism_json(data=data,
                                                     file_uuid=generate_file_uuid(bundle_uuid, file_name),
                                                     i=donor_number)
    with open(f'{output_dir}/{file_name}', 'w') as f:
        f.write(json.dumps(donor_organism_json, indent=4))
    print(f'"{output_dir}/{file_name}" successfully written.')


def run(xlsx, output_dir, upload=False):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    wb = load_workbook(xlsx)

    project_json, project_uuid = generate_project_json(wb, output_dir)
    bundle_uuid = copy.deepcopy(project_uuid)

    cell_counts = get_cell_counts()
    cell_count = 1  # the default if not found
    for accession in project_json['geo_series_accessions']:
        if accession in cell_counts:
            cell_count = cell_counts[accession]


    # add_matrix_file(project_json['geo_series_accessions'], project_uuid, output_dir)
    generate_cell_suspension_json(wb=wb,
                                  output_dir=output_dir,
                                  cell_count=cell_count,
                                  bundle_uuid=bundle_uuid)
    generate_specimen_from_organism_json(wb=wb,
                                         output_dir=output_dir,
                                         bundle_uuid=bundle_uuid)
    generate_donor_organism_jsons(wb=wb,
                                  output_dir=output_dir,
                                  bundle_uuid=bundle_uuid)
    generate_links_json(output_dir)

    if upload:
        hca_config = HCAConfig()

        hca_config["DSSClient"].swagger_url = f"https://dss.dev.data.humancellatlas.org/v1/swagger.json"
        dss = DSSClient(config=hca_config)

        response = dss.upload(src_dir=output_dir,
                              replica='aws',
                              staging_bucket='lon-test-data',
                              bundle_uuid=bundle_uuid)
        print(f'Successful upload.  Bundle information is:\n{json.dumps(response, indent=4)}')


def main(argv):
    parser = argparse.ArgumentParser(description='Turn an xlsx file into the json files necessary for a complete '
                                                 'bundle to be uploaded into the DSS and serve as a minimal project.')
    parser.add_argument("--xlsx", type=str,
                        help="Path to an xlsx (excel) file.  "
                             "Example: 'data/test_project_000.xlsx'")
    parser.add_argument("--output_dir", type=str,
                        default='bundle',
                        help="Path to an output directory.")
    parser.add_argument("--upload", type=bool,
                        default=False,
                        help="Whether or not one should upload this data as a bundle to the data-store.")

    args = parser.parse_args(argv)
    run(args.xlsx, args.output_dir, args.upload)


if __name__ == "__main__":
    main(sys.argv[1:])
