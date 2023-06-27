import subprocess
import json
import os


def run_scaffold(main, file_path, output_path):
    command = f"python kmer_main.py -f {main} -p {file_path} -o {output_path}"
    print(f"Running command: {command}")
    subprocess.call(command, shell=True)
    print(f"Command completed: {command}")

def run_kmer_extraction(main, file_path, output_path, k_value, foldername, core_num):
    command = f"python kmer_main.py -f {main} -k {k_value} -p {file_path} -o {output_path} -fn {foldername} -c {core_num}"
    print(f"Running command: {command}")
    subprocess.call(command, shell=True)
    print(f"Command completed: {command}")


json_file = "start_data copy.json"
json_path = os.path.join(os.getcwd(), json_file)
json_data = open(json_path)



data = json.load(json_data)

# Variables to store the values
k_mer = None
path = None
output = None
func = None
foldername = None
cores = None

# Dictionaries to store the values
path_list = {}
foldername_list = {}
func_list = {}

# Assign values to variables or dictionaries based on the context
for key, value in data['scaffold'].items():
    if isinstance(value['value'], list):
        if key == 'path':
            path_list = value['value']
        elif key =='func':
            func_list = value['value']
    else:
        if key == 'output':
            output = value['value']



total_commands = len(path_list) 
completed_commands = 0

for i, func_name in enumerate(func_list):
    path = path_list[i]
    run_scaffold(func_name, path, output)
    completed_commands += 1
    progress = (completed_commands / total_commands) * 100
    print(f"Progress: {progress:.2f}%")



# Assign values to variables or dictionaries based on the context
for key, value in data['kmer'].items():
    if isinstance(value['value'], list):
        if key == 'path':
            path_list = value['value']
        elif key == 'foldername':
            foldername_list = value['value']
    else:
        if key == 'k_mer':
            k_mer = value['value']
        elif key == 'output':
            output = value['value']
        elif key == 'func':
            func = value['value']
        elif key == 'cores':
            cores = int(value['value'])


total_commands = len(path_list) 
completed_commands = 0

for i, path_name in enumerate(path_list):
    foldername = foldername_list[i]
    run_kmer_extraction(func, path_name, output, k_mer, foldername, cores)
    completed_commands += 1
    progress = (completed_commands / total_commands) * 100
    print(f"Progress: {progress:.2f}%")


