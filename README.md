# load-project

This will take an xlsx file and generate a `project_0.json` suitable for uploading to the DSS.

This will then upload the `project_0.json` and a (mostly) empty `links.json` which will populate a new project in 
the browser.

Source a new environment and install dependencies:

    virtualenv -p python3.6 v3nv && . v3nv/bin/activate && pip install -r requirements.txt

Parse the xlsx:

    #!/usr/bin/env bash

    # E-GEOD-81547_curated_ontologies_07_2019.xlsx
    python xlsx_to_project_json.py --xlsx data/test_000.xlsx
    
    # Gary_Bader_9_16.xlsx
    python xlsx_to_project_json.py --xlsx data/test_001.xlsx
    
    # GEOD-93593_HCA_Ontologies_July_2.xlsx
    python xlsx_to_project_json.py --xlsx data/test_002.xlsx
    
    # hca-metadata-spreadsheet-GSE84133_pancreas.xlsx
    python xlsx_to_project_json.py --xlsx data/test_003.xlsx
    
    # hca-metadata-spreadsheet-GSE95459-GSE114374-colon.xlsx
    python xlsx_to_project_json.py --xlsx data/test_004.xlsx
    
    # mf-E-GEOD-106540_spreadsheet_v9.xlsx
    python xlsx_to_project_json.py --xlsx data/test_005.xlsx

Adding `--upload true` will upload the data to the DSS.
Note that UUID's are now always programmatically generated from GEO accessions 
and cannot be provided via the commandline.

**NOTE**:

Edited the following fields in "data/test_004.xlsx":

    publications.publication_url -> publications.url
    publications.publication_title -> publications.title
