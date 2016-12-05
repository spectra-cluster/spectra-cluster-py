all: clean dist

init:
	pip3 install -r requirements.txt

clean:
	rm -rf build dist

test:
	cd tests && \
		python3 -m unittest

dist:
	pyinstaller --clean --onefile spectra_cluster/ui/protein_annotator.py
	pyinstaller --clean --onefile spectra_cluster/ui/id_transferer_cli.py
	pyinstaller --clean --onefile spectra_cluster/ui/mgf_search_result_annotator.py
	pyinstaller --clean --onefile spectra_cluster/ui/cluster_features_cli.py
	pyinstaller --clean --onefile spectra_cluster/ui/consensus_spectrum_exporter.py

.PHONY:
	clean dist
