# load-project

This will take an xlsx file and generate a `project_0.json` suitable for uploading to the DSS.

This will then upload the `project_0.json` and a (mostly) empty `links.json` which will populate a new project in 
the browser.

Source a new environment and install dependencies:

    virtualenv -p python3.6 v3nv && . v3nv/bin/activate && pip install -r requirements.txt

Parse the xlsx:

    #!/usr/bin/env bash

    # E-GEOD-81547_curated_ontologies_07_2019.xlsx
    # DSS prod uuid: cddab57b-6868-4be4-806f-395ed9dd635a
    python xlsx_to_project_json.py --xlsx data/test_000.xlsx
    
    # Gary_Bader_9_16.xlsx
    # DSS prod uuid: 4d6f6c96-2a83-43d8-8fe1-0f53bffd4674
    python xlsx_to_project_json.py --xlsx data/test_001.xlsx
    
    # GEOD-93593_HCA_Ontologies_July_2.xlsx
    # DSS prod uuid: 2043c65a-1cf8-4828-a656-9e247d4e64f1
    python xlsx_to_project_json.py --xlsx data/test_002.xlsx
    
    # hca-metadata-spreadsheet-GSE84133_pancreas.xlsx
    # DSS prod uuid: f86f1ab4-1fbb-4510-ae35-3ffd752d4dfc
    python xlsx_to_project_json.py --xlsx data/test_003.xlsx
    
    # hca-metadata-spreadsheet-GSE95459-GSE114374-colon.xlsx
    # DSS prod uuid: f8aa201c-4ff1-45a4-890e-840d63459ca2
    python xlsx_to_project_json.py --xlsx data/test_004.xlsx
    
    # mf-E-GEOD-106540_spreadsheet_v9.xlsx
    # DSS prod uuid: 90bd6933-40c0-48d4-8d76-778c103bf545
    python xlsx_to_project_json.py --xlsx data/test_005.xlsx

Adding `--upload true` will upload the data to the DSS.
Note that UUID's are now always programmatically generated from GEO accessions 
and cannot be provided via the commandline.

**NOTE**:

Edited the following fields in "data/test_004.xlsx":

    publications.publication_url -> publications.url
    publications.publication_title -> publications.title
