import argparse
import copy
from datetime import datetime
import json
import os
from pathlib import Path
import sys
from typing import (
    List,
    Sequence,
)
import uuid

from openpyxl import load_workbook

from count_cells import get_cell_counts


# TODO: Consolidate similar functions and clean up code.


def timestamp():
    return datetime.utcnow().strftime("%Y-%m-%dT%H%M%S.%fZ")


def generate_project_uuid(geo_accessions: Sequence[str]) -> str:
    """
    Deterministically generate a project UUID based on one or more GEO accession ids.
    """
    if isinstance(geo_accessions, str):
        geo_accessions = [geo_accessions]
    namespace_uuid = uuid.UUID('296e1002-1c99-4877-bb7f-bb6a3b752638')
    return str(uuid.uuid5(namespace_uuid, ''.join(sorted(geo_accessions))))


def generate_file_uuid(bundle_uuid: str, file_name: str) -> str:
    """
    Deterministically generate a file UUID based on the parent bundle uuid and its file name.
    """
    namespace_uuid = uuid.UUID('4c52e3d0-ffe5-4b4d-a4d0-cb6a6f372b31')
    return str(uuid.uuid5(namespace_uuid, bundle_uuid + file_name))


def get_spreadsheet_paths(src_dir: Path) -> List[Path]:
    return [p for p in src_dir.iterdir() if p.is_file() and p.name.endswith('.0.xlsx')]


def create_project_json(data, version, verify=False):
    project_json = {
      "describedBy": "https://schema.dev.data.humancellatlas.org/type/project/14.0.0/project",
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
    funders = []
    for i in range(len(grant_ids)):
        funder = {}
        if data["funders.grant_title"][i]:
            funder['grant_title'] = str(data["funders.grant_title"][i])
        if data["funders.grant_id"][i]:
            funder['grant_id'] = str(data["funders.grant_id"][i])
        if data["funders.organization"][i]:
            funder['organization'] = str(data["funders.organization"][i])
        if funder:
            funders.append(funder)
    if funders:
        project_json['funders'] = funders
    else:
        project_json['funders'] = [
            {
                'grant_title': 'none given',
                'grant_id': 'none given',
                'organization': 'none given'
            }
        ]

    project_uuid = generate_project_uuid(project_json['geo_series_accessions'])
    project_json["provenance"] = {
        "document_id": project_uuid,
        "submission_date": version,
        "update_date": version
        }
    if verify:
        print(json.dumps(project_json, indent=4))
    return project_json, project_uuid


def create_cell_suspension_jsons(data, cell_count, file_uuid, i=0):
    version = timestamp()
    cell_suspension_json = {
        "describedBy": "https://schema.dev.data.humancellatlas.org/type/biomaterial/13.1.0/cell_suspension",
        "schema_type": "biomaterial",
        'biomaterial_core':
            {
                'biomaterial_id': None,
                'biomaterial_name': None,
                'biomaterial_description': None,
                'ncbi_taxon_id': [0],
                'genotype': None
            },
        'genus_species':
            [
                {
                    'text': None,
                    'ontology': None,
                    'ontology_label': None
                }
            ]
    }

    cell_suspension_json.update(
        fill_sections(
            data=data,
            as_list=[
                'genus_species',
                'selected_cell_types'],
            keys={
                'biomaterial_core':
                    [
                        'biomaterial_id',
                        'biomaterial_name',
                        'biomaterial_description',
                        'ncbi_taxon_id',
                        # 'biosamples_accession',
                        # 'insdc_sample_accession',
                        'genotype'
                    ],
                'genus_species':
                    [
                        'text',
                        'ontology',
                        'ontology_label'
                    ],
                # 'cell_morphology':
                #     [
                #         'percent_cell_viability',
                #         'cell_viability_method'
                #     ],
                'selected_cell_types':
                    [
                        'text',
                        'ontology',
                        'ontology_label'
                    ],
                'timecourse':
                    [
                        'value',
                        'unit'
                    ],
            },
            i=i
        )
    )
    cell_suspension_json['estimated_cell_count'] = cell_count
    timecourse_unit = conditionally_create_subdict(
        data=data,
        primary_key='timecourse.unit',
        secondary_keys=[
            'text',
            'ontology',
            'ontology_label'
        ],
        i=i,
        use_list=False
    )
    if timecourse_unit:
        cell_suspension_json.get('timecourse', {}).update(timecourse_unit)

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
                    if primary_key == 'paired_end':
                        if str(data[primary_key][i]) == 'yes':
                            document[primary_key] = True
                        elif str(data[primary_key][i]) == 'no':
                            document[primary_key] = False
                        else:
                            raise RuntimeError('This should never happen.')
                    else:
                        document[primary_key] = str(data[primary_key][i])
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
                    sub_dictionary[key] = str(data[f'{primary_key}.{key}'][i])
    if sub_dictionary:
        prime_dictionary[primary_key.split('.')[-1]] = [sub_dictionary] if use_list else sub_dictionary
    return prime_dictionary


def create_sequencing_protocol_json(data, file_uuid, i=0):
    version = timestamp()
    sequencing_protocol_json = {
        "describedBy": "https://schema.dev.data.humancellatlas.org/type/protocol/sequencing/10.0.0/sequencing_protocol",
        "schema_type": "protocol",
        'protocol_core':
            {
                'protocol_id': None,
                'protocol_name': None,
                'protocol_description': None
            },
        'instrument_manufacturer_model':
            {
                'text': None,
                'ontology': None,
                'ontology_label': None
            },
        'paired_end': None,
        'method':
            {
                'text': None,
                'ontology': None,
                'ontology_label': None
            }
    }
    sequencing_protocol_json.update(
        fill_sections(
            data=data,
            as_list=[],
            keys={
                'protocol_core':
                    [
                        'protocol_id',
                        'protocol_name',
                        'protocol_description'
                    ],
                'instrument_manufacturer_model':
                    [
                        'text',
                        'ontology',
                        'ontology_label'
                    ],
                'paired_end': [],
                'method':
                    [
                        'text',
                        'ontology',
                        'ontology_label'
                    ]
            },
            i=i
        )
    )
    sequencing_protocol_json["provenance"] = {
        "document_id": file_uuid,
        "submission_date": version,  # TODO: Fetch from DSS if it exists
        "update_date": version
    }
    return sequencing_protocol_json


def create_library_preparation_protocol_json(data, file_uuid, i=0):
    version = timestamp()
    library_preparation_protocol_json = {
        "describedBy": "https://schema.dev.data.humancellatlas.org/type/protocol/sequencing/6.1.0/library_preparation_protocol",
        "schema_type": "protocol",
        'protocol_core':
            {
                'protocol_id': None,
                'protocol_name': None,
                'protocol_description': None
            },
        'library_construction_method':
            {
                'text': None,
                'ontology': None,
                'ontology_label': None
            },
        'input_nucleic_acid_molecule':
            {
                'text': None,
                'ontology': None,
                'ontology_label': None
            },
        'nucleic_acid_source': None,
        'end_bias': None,
        'strand': None
    }
    library_preparation_protocol_json.update(
        fill_sections(
            data=data,
            as_list=[],
            keys={
                'protocol_core':
                    [
                        'protocol_id',
                        'protocol_name',
                        'protocol_description'
                    ],
                'nucleic_acid_source': [],
                'input_nucleic_acid_molecule':
                    [
                        'text',
                        'ontology',
                        'ontology_label'
                    ],
                'library_construction_method':
                    [
                        'text',
                        'ontology',
                        'ontology_label'
                    ],
                'library_construction_kit':
                    [
                        'retail_name',
                        'catalog_number',
                        'manufacturer'
                    ],
                'end_bias': [],
                'primer': [],
                'strand': []
            },
            i=i
        )
    )
    library_preparation_protocol_json["provenance"] = {
        "document_id": file_uuid,
        "submission_date": version,  # TODO: Fetch from DSS if it exists
        "update_date": version
    }
    return library_preparation_protocol_json


def create_specimen_from_organism_json(data, file_uuid, i=0):
    version = timestamp()
    specimen_from_organism_json = {
        "describedBy": "https://schema.dev.data.humancellatlas.org/type/biomaterial/10.2.0/specimen_from_organism",
        "schema_type": "biomaterial",
        'biomaterial_core':
            {
                'biomaterial_id': None,
                'biomaterial_name': None,
                'biomaterial_description': None,
                'ncbi_taxon_id': [0],
                'genotype': None
            },
        'organ':
            {
                'text': None,
                'ontology': None,
                'ontology_label': None
            },
        'genus_species':
            [
                {
                    'text': None,
                    'ontology': None,
                    'ontology_label': None
                }
            ]
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
                        # 'biosamples_accession',
                        # 'insdc_sample_accession',
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
        raise RuntimeError(f'This should never happen.  ({type(ids)})')


def create_donor_organism_json(data, file_uuid, i=0):
    version = timestamp()
    donor_organism_json = {
        "describedBy": "https://schema.dev.data.humancellatlas.org/type/biomaterial/15.3.0/donor_organism",
        "schema_type": "biomaterial",
        'biomaterial_core':
            {
                'biomaterial_id': None,
                'biomaterial_name': None,
                'biomaterial_description': None,
                'ncbi_taxon_id': [0],
                'genotype': None
            },
        'development_stage':
            {
                'text': None,
                'ontology': None,
                'ontology_label': None
            },
        'genus_species':
            [
                {
                    'text': None,
                    'ontology': None,
                    'ontology_label': None
                }
            ],
        'is_living': None,
        'sex': None
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
                        # 'biosamples_accession',
                        # 'insdc_sample_accession',
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
                # Weird parsing???
                # https://logs.dev.data.humancellatlas.org/_plugin/kibana/app/kibana#/discover?_g=(refreshInterval:(display:Off,pause:!f,value:0),time:(from:'2019-12-20T05:51:59.104Z',mode:absolute,to:'2019-12-20T05:51:59.456Z'))&_a=(columns:!(_source),filters:!(('$state':(store:appState),meta:(alias:!n,disabled:!f,index:'*',key:'@log_group',negate:!f,type:phrase,value:%2Faws%2Flambda%2Fdss-index-dev),query:(match:('@log_group':(query:%2Faws%2Flambda%2Fdss-index-dev,type:phrase))))),index:'*',interval:auto,query:(query_string:(analyze_wildcard:!t,query:ERROR)),sort:!('@timestamp',desc))
                # 'organism_age': [],
                # 'organism_age_unit':
                #     [
                #         'text',
                #         'ontology',
                #         'ontology_label'
                #     ],
                # 'human_specific':
                #     [
                #         'body_mass_index'
                #     ],
                # 'death':
                #     [
                #         'cause_of_death'
                #     ],
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
    # ethnicity = conditionally_create_subdict(
    #     data=data,
    #     primary_key='human_specific.ethnicity',
    #     secondary_keys=[
    #         'text',
    #         'ontology',
    #         'ontology_label'
    #     ],
    #     i=i,
    #     use_list=True
    # )
    # if ethnicity:
    #     donor_organism_json.get('human_specific', {}).update(ethnicity)
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


def parse_library_preparation_protocol_data_from_xlsx(wb):
    data = {}

    section_keywords = ['library']
    project_sections = get_harmonized_project_sections(wb, section_keywords)
    for project_key in section_keywords:
        try:
            data.update(process_section(section=project_sections[project_key], prefix='library_preparation_protocol.'))
        except KeyError:
            print(f'{project_key} section not found in the data!')
    return data


def parse_sequencing_protocol_data_from_xlsx(wb):
    data = {}

    section_keywords = ['sequencing']
    project_sections = get_harmonized_project_sections(wb, section_keywords)
    for project_key in section_keywords:
        try:
            data.update(process_section(section=project_sections[project_key], prefix='sequencing_protocol.'))
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


def write_project_json(project_json, output_dir):
    with open(f'{output_dir}/project_0.json', 'w') as f:
        f.write(json.dumps(project_json, indent=4))
    print(f'"{output_dir}/project_0.json" successfully written.')


def generate_cell_suspension_json(wb, output_dir, cell_count, bundle_uuid):
    file_name = 'cell_suspension_0.json'
    cell_suspension_data = parse_cell_suspension_data_from_xlsx(wb)
    cell_json = create_cell_suspension_jsons(data=cell_suspension_data,
                                             cell_count=cell_count,
                                             file_uuid=generate_file_uuid(bundle_uuid, file_name))
    with open(f'{output_dir}/{file_name}', 'w') as f:
        f.write(json.dumps(cell_json, indent=4))
    print(f'"{output_dir}/{file_name}" successfully written.')


def generate_specimen_from_organism_jsons(wb, output_dir, bundle_uuid):
    data = parse_specimen_from_organism_data_from_xlsx(wb)
    # specimens = [specimen for specimen in data['biomaterial_core.biomaterial_id'] if specimen]
    file_number = 0
    seen_organs = set()
    for specimen_number, organ in enumerate(data['organ.text']):
        if organ not in seen_organs:
            if organ:
                generate_specimen_from_organism_json(data, output_dir, specimen_number, file_number, bundle_uuid)
                file_number += 1
            seen_organs.add(organ)
    if file_number == 0:
        generate_specimen_from_organism_json(data, output_dir, 0, 0, bundle_uuid)


def generate_specimen_from_organism_json(data, output_dir, specimen_number, file_number, bundle_uuid):
    file_name = f'specimen_from_organism_{file_number}.json'
    specimen_from_organism_json = create_specimen_from_organism_json(
        data=data,
        file_uuid=generate_file_uuid(bundle_uuid, file_name),
        i=specimen_number)
    with open(f'{output_dir}/{file_name}', 'w') as f:
        f.write(json.dumps(specimen_from_organism_json, indent=4))
    print(f'"{output_dir}/{file_name}" successfully written.')


def generate_library_preparation_protocol_json(wb, output_dir, bundle_uuid):
    file_name = 'library_preparation_protocol_0.json'
    library_preparation_protocol_data = parse_library_preparation_protocol_data_from_xlsx(wb)
    library_preparation_protocol_json = create_library_preparation_protocol_json(
        data=library_preparation_protocol_data,
        file_uuid=generate_file_uuid(bundle_uuid, file_name)
    )
    with open(f'{output_dir}/{file_name}', 'w') as f:
        f.write(json.dumps(library_preparation_protocol_json, indent=4))
    print(f'"{output_dir}/{file_name}" successfully written.')


def generate_analysis_protocol_json(output_dir, bundle_uuid):
    # TODO: Hard-coded and not sure where this data should come from... ???
    file_name = 'analysis_protocol_0.json'
    version = timestamp()
    analysis_protocol_json = {
        "computational_method": "SmartSeq2SingleCell",
        "describedBy": "https://schema.humancellatlas.org/type/protocol/analysis/9.0.0/analysis_protocol",
        "protocol_core": {
            "protocol_id": "smartseq2_v2.3.0"
        },
        "schema_type": "protocol",
        "type": {
            "text": "analysis"
        },
        "provenance": {
            "document_id": generate_file_uuid(bundle_uuid, file_name),
            "submission_date": version,  # TODO: Fetch from DSS if it exists
            "update_date": version
        }
    }
    with open(f'{output_dir}/{file_name}', 'w') as f:
        f.write(json.dumps(analysis_protocol_json, indent=4))
    print(f'"{output_dir}/{file_name}" successfully written.')


def generate_sequencing_protocol_json(wb, output_dir, bundle_uuid):
    file_name = 'sequencing_protocol_0.json'
    sequencing_protocol_data = parse_sequencing_protocol_data_from_xlsx(wb)
    sequencing_protocol_json = create_sequencing_protocol_json(
        data=sequencing_protocol_data,
        file_uuid=generate_file_uuid(bundle_uuid, file_name)
    )
    with open(f'{output_dir}/{file_name}', 'w') as f:
        f.write(json.dumps(sequencing_protocol_json, indent=4))
    print(f'"{output_dir}/{file_name}" successfully written.')


def generate_analysis_json(bundle_uuid, output_dir):
    file_name = 'analysis_file_0.json'
    version = timestamp()
    analysis_json = {
        "describedBy": "https://schema.humancellatlas.org/type/file/6.0.0/analysis_file",
        "file_core": {
            "file_name": "matrix.mtx.zip",
            "format": "mtx"
        },
        "schema_type": "file",
        "provenance": {
            "document_id": generate_file_uuid(bundle_uuid, file_name),
            "submission_date": version,  # TODO: Fetch from DSS if it exists
            "update_date": version
        }
    }
    with open(f'{output_dir}/{file_name}', 'w') as f:
        f.write(json.dumps(analysis_json, indent=4))
    print(f'"{output_dir}/{file_name}" successfully written.')


def generate_links_json(output_dir):
    file_name = 'links.json'
    links_json = {
        'describedBy': 'https://schema.humancellatlas.org/system/1.1.5/links',
        'schema_type': 'link_bundle',
        'schema_version': '1.1.5',
        'links': []
    }
    with open(f'{output_dir}/{file_name}', 'w') as f:
        f.write(json.dumps(links_json, indent=4))
    print(f'"{output_dir}/{file_name}" successfully written.')


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


def run(xlsx, output_dir=None, clear=True):
    wb = load_workbook(xlsx)

    project_data = parse_project_data_from_xlsx(wb)
    project_json, project_uuid = create_project_json(project_data, version=timestamp())

    root = f'projects/{project_uuid}'
    matrix_file = f'{root}/bundle/matrix.mtx.zip'
    output_dir = f'{root}/bundle' if not output_dir else output_dir

    if clear and os.path.exists(output_dir):
        remove_previous_metadata(output_dir=output_dir)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    write_project_json(project_json, output_dir)
    bundle_uuid = copy.deepcopy(project_uuid)

    if os.path.exists(matrix_file):
        generate_analysis_json(bundle_uuid=bundle_uuid, output_dir=output_dir)

    cell_counts = get_cell_counts()
    cell_count = 1  # the default if not found
    for accession in project_json['geo_series_accessions']:
        if accession in cell_counts:
            cell_count = cell_counts[accession]

    generate_cell_suspension_json(wb=wb,
                                  output_dir=output_dir,
                                  cell_count=cell_count,
                                  bundle_uuid=bundle_uuid)
    generate_specimen_from_organism_jsons(wb=wb,
                                          output_dir=output_dir,
                                          bundle_uuid=bundle_uuid)
    generate_donor_organism_jsons(wb=wb,
                                  output_dir=output_dir,
                                  bundle_uuid=bundle_uuid)
    generate_library_preparation_protocol_json(wb=wb,
                                               output_dir=output_dir,
                                               bundle_uuid=bundle_uuid)
    generate_sequencing_protocol_json(wb=wb,
                                      output_dir=output_dir,
                                      bundle_uuid=bundle_uuid)
    # generate_analysis_protocol_json(output_dir=output_dir,
    #                                 bundle_uuid=bundle_uuid)
    generate_links_json(output_dir)


def remove_previous_metadata(output_dir):
    for root, dirs, files in os.walk(output_dir):
        for name in files:
            if name.strip().endswith('.json'):
                os.remove(os.path.join(root, name))
    print(f'All previous jsons cleared from {output_dir}.')


def main(argv):
    parser = argparse.ArgumentParser(description='Turn an xlsx file into the json files necessary for a complete '
                                                 'bundle to be uploaded into the DSS and serve as a minimal project.')
    parser.add_argument("--xlsx", type=str,
                        help="Path to an xlsx (excel) file.  "
                             "Example: 'data/test_project_000.xlsx'")
    parser.add_argument("--output_dir", type=str,
                        help="Path to an output directory.")
    parser.add_argument("--no-clean",
                        dest='clean', action='store_false', required=False,
                        help="Whether not to delete already present jsons from the output directory "
                             "prior to populating.")
    parser.add_argument("--clean",
                        dest='clean', action='store_true', required=False,
                        help="Whether to delete already present jsons from the output directory prior to populating.")
    parser.set_defaults(clean=True)

    args = parser.parse_args(argv)

    run(args.xlsx, args.output_dir, args.clean)


if __name__ == "__main__":
    main(sys.argv[1:])
