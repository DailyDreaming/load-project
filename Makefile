all:
	$(MAKE) download
	$(MAKE) extract
	$(MAKE) merge
	$(MAKE) count
	$(MAKE) metadata
	$(MAKE) bundles
	$(MAKE) assets

download:
	python3 download.py --all

extract:
	python3 extract.py

merge:
	@echo 'Replace with merge matrix script'

count:
	python3 count_cells.py --write-all

metadata:
	python3 generate_metadata.py

bundles:
	@echo 'Replace with upload bundles script'

assets:
	@echo 'Replace with upload to assets bucket script'

.PHONY: all download extract merge count metadata bundles assets
