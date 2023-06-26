import subprocess
import json
import os


def run_kmer_extraction(main, file_path, output_path, k_value, foldername, core_num):
    command = f"python kmer_main.py -f {main} -k {k_value} -p {file_path} -o {output_path} -fn {foldername} -c {core_num}"
    print(f"Running command: {command}")
    subprocess.call(command, shell=True)
    print(f"Command completed: {command}")



json_file = "start_data.json"
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
for key, value in data.items():
    if isinstance(value['value'], list):
        if key == 'path':
            path_list = value['value']
        elif key == 'foldername':
            foldername_list = value['value']
        elif key =='func':
            func_list = value['value']
    else:
        if key == 'k_mer':
            k_mer = value['value']
        elif key == 'output':
            output = value['value']
        elif key == 'cores':
            cores = int(value['value'])

# Print the assigned values
print("k_mer:", k_mer)
print("output:", output)
print("cores:", cores)
# Print the path dictionary
print("Path List:", str(path_list))
# Print the foldername dictionary
print("Foldername List:", str(foldername_list))
# Print the func dictionary
print("Func List:", str(func_list))


total_commands = len(path_list) 
completed_commands = 0

for i, file_path in enumerate(path_list):
    foldername = foldername_list[i]
    if foldername =='wgs':
        run_kmer_extraction(func_list[0], file_path, output, k_mer, foldername, cores)
    else:
        run_kmer_extraction(func_list[1],file_path, output, k_mer, foldername, cores)
    completed_commands += 1
    progress = (completed_commands / total_commands) * 100
    print(f"Progress: {progress:.2f}%")
