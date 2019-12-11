#!/usr/bin/env bash

virtualenv -p python3.6 v3nv && . v3nv/bin/activate && pip install -r requirements.txt && python xlsx_to_project_json.py --uuid 4d6f6c96-2a83-43d8-8fe1-0f53bffd4674 --xlsx Gary_Bader_9_16.xlsx
cat project_0.json
