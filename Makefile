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
	python script_boost.py -j start_data_boost.json
	python plot.py


scaffold:
	python script_scaffold.py -j start_data_scaffold.json

kmer:
	python script_kmer.py -j start_data_kmer.json

gexp:
	python script_gexp.py -j start_data_gexp.json

target:
	python script_target.py -j start_data_target.json

models:
	python script_models.py -j start_data_models.json

boost:
	python script_boost.py -j start_data_boost.json

plot:
	python plot.py


scaffoldunix:
	python script_scaffold.py -j start_data_models_unix.json

kmerunix:
	python script_kmer.py -j start_data_models_unix.json

gexpunix:
	python script_gexp.py -j start_data_models_unix.json

targetunix:
	python script_target.py -j start_data_target_unix.json

modelsunix:
	python script_models.py -j start_data_models_unix.json

boostunix:
	python script_boost.py -j start_data_boost_unix.json
