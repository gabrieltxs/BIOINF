# Makefile

.PHONY: install run scaffold kmer gexp target modelswin modelsunix

install:
	pip install -r requirements.txt

run:
	python script_scaffold.py -j start_data_scaffold.json
	python script_kmer.py -j start_data_kmer.json
	python script_gexp.py -j start_data_gexp.json
	python script_target.py -j start_data_target.json
	python script_models.py -j start_data_models.json


scaffold:
	python script_scaffold.py -j start_data_scaffold.json

kmer:
	python script_kmer.py -j start_data_kmer.json

gexp:
	python script_gexp.py -j start_data_gexp.json

target:
	python script_target.py -j start_data_target.json

modelswin:
	python script_models.py -j start_data_models.json

modelsunix:
	python script_models.py -j start_data_models_unix.json
