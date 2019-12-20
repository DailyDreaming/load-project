all:
	$(MAKE) download
	$(MAKE) extract
	$(MAKE) merge
	$(MAKE) count
	$(MAKE) metadata
	$(MAKE) bundles
	$(MAKE) assets

download:
	python3 ./download_geo_matrix.py

extract:
	python3 ./unpack_tar_file.py

merge:
	@echo 'Replace with merge matrix script'

count:
	@echo 'Replace with merge count cells script'

metadata:
	python3 ./run_all.py

bundles:
	@echo 'Replace with upload bundles script'

assets:
	@echo 'Replace with upload to assets bucket script'

.PHONY: all download extract merge count metadata bundles assets
