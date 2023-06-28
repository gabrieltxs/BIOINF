import subprocess
import json
import os




def run_gexp(filename, list_genomes):
    command = f"python gexp.py -f {filename} -l {list_genomes}"
    print(f"Running command: {command}")
    subprocess.call(command, shell=True)
    print(f"Command completed: {command}")


json_file = "start_data_gexp.json"
json_path = os.path.join(os.getcwd(), json_file)
json_data = open(json_path)



data = json.load(json_data)



# Assign values to variables or dictionaries based on the context
for key, value in data.items():
    if key == 'filename':
        filename = value['value']
    elif key == 'list_genomes':
        list_genomes = value['value']

run_gexp(filename, list_genomes)
    


