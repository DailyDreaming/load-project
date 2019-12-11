# load-project

This will take an xlsx file and generate a `project_o.json` suitable for uploading to the DSS.

This will then upload the `project_o.json` and a (mostly) empty `link.json` which will populate a new project in 
the browser.

Source a new environment and install dependencies:

    virtualenv -p python3.6 v3nv && . v3nv/bin/activate && pip install -r requirements.txt

Parse the xlsx and include a project uuid:

    python xlsx_to_project_json.py --uuid 4d6f6c96-2a83-43d8-8fe1-0f53bffd4674 --xlsx Gary_Bader_9_16.xlsx

Upload to the DSS:

    python upload_to_dss.py
