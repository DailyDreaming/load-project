#!/usr/bin/env bash

# E-GEOD-81547_curated_ontologies_07_2019.xlsx
python xlsx_to_project_json.py --uuid cddab57b-6868-4be4-806f-395ed9dd635a --xlsx data/test_000.xlsx
python upload_to_dss.py

# Gary_Bader_9_16.xlsx
python xlsx_to_project_json.py --uuid 4d6f6c96-2a83-43d8-8fe1-0f53bffd4674 --xlsx data/test_001.xlsx
python upload_to_dss.py

# GEOD-93593_HCA_Ontologies_July_2.xlsx
python xlsx_to_project_json.py --uuid 2043c65a-1cf8-4828-a656-9e247d4e64f1 --xlsx data/test_002.xlsx
python upload_to_dss.py

# hca-metadata-spreadsheet-GSE84133_pancreas.xlsx
python xlsx_to_project_json.py --uuid f86f1ab4-1fbb-4510-ae35-3ffd752d4dfc --xlsx data/test_003.xlsx
python upload_to_dss.py

# hca-metadata-spreadsheet-GSE95459-GSE114374-colon.xlsx
python xlsx_to_project_json.py --uuid f8aa201c-4ff1-45a4-890e-840d63459ca2 --xlsx data/test_004.xlsx
python upload_to_dss.py

# mf-E-GEOD-106540_spreadsheet_v9.xlsx
python xlsx_to_project_json.py --uuid 90bd6933-40c0-48d4-8d76-778c103bf545 --xlsx data/test_005.xlsx
python upload_to_dss.py