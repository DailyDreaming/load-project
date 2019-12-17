#!/usr/bin/env bash

# E-GEOD-81547_curated_ontologies_07_2019.xlsx
python xlsx_to_project_json.py --xlsx data/test_000.xlsx
python upload_to_dss.py

# Gary_Bader_9_16.xlsx
python xlsx_to_project_json.py --xlsx data/test_001.xlsx
python upload_to_dss.py

# GEOD-93593_HCA_Ontologies_July_2.xlsx
python xlsx_to_project_json.py --xlsx data/test_002.xlsx
python upload_to_dss.py

# hca-metadata-spreadsheet-GSE84133_pancreas.xlsx
python xlsx_to_project_json.py --xlsx data/test_003.xlsx
python upload_to_dss.py

# hca-metadata-spreadsheet-GSE95459-GSE114374-colon.xlsx
python xlsx_to_project_json.py --xlsx data/test_004.xlsx
python upload_to_dss.py

# mf-E-GEOD-106540_spreadsheet_v9.xlsx
python xlsx_to_project_json.py --xlsx data/test_005.xlsx
python upload_to_dss.py
