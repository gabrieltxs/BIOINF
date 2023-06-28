# Makefile

.PHONY: install run

install:
	pip install -r requirements.txt

run:
	python script_scaffold.py -j start_data_scaffold.json
