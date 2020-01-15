all:
	$(MAKE) download
	$(MAKE) extract
	$(MAKE) merge
	$(MAKE) count
	$(MAKE) metadata
	$(MAKE) bundles
	$(MAKE) assets

download:
	python3 download.py

extract:
	python3 extract.py

matrices:
	python3 convert_matrices.py

clean_matrices:
	python clean.py

count:
	python3 count_cells.py --write

metadata:
	python3 generate_metadata.py

bundles:
	python3 upload_bundles.py

assets:
	python3 upload_assets.py

report:
	python3 overview_report.py | column -t | less -S

.PHONY: all download extract merge count metadata bundles assets
