import subprocess
import json
import os
import argparse


def run_kmer_extraction(main, file_path, output_path, k_value, foldername, core_num):
    command = f"python kmer.py -f {main} -k {k_value} -p {file_path} -o {output_path} -fn {foldername} -c {core_num}"
    print(f"Running command: {command}")
    subprocess.call(command, shell=True)
    print(f"Command completed: {command}")


def main(json_file):
    json_path = os.path.join(os.getcwd(), json_file)
    with open(json_path) as json_data:
        data = json.load(json_data)

    # Variables to store the values
    k_mer = None
    output = None
    cores = None

    # Lists to store the values
    path_list = []
    foldername_list = []
    func_list = []

    # Assign values to variables or lists based on the context
    for key, value in data.items():
        if isinstance(value['value'], list):
            if key == 'path':
                path_list = value['value']
            elif key == 'foldername':
                foldername_list = value['value']
            elif key == 'func':
                func_list = value['value']
        else:
            if key == 'k_mer':
                k_mer = value['value']
            elif key == 'output':
                output = value['value']
            elif key == 'cores':
                cores = int(value['value'])

    total_commands = len(path_list)
    completed_commands = 0

    for i, path_name in enumerate(path_list):
        foldername = foldername_list[i]
        if foldername == 'wgs':
            run_kmer_extraction(func_list[1], path_name, output, k_mer, foldername, cores)
        else:
            run_kmer_extraction(func_list[0], path_name, output, k_mer, foldername, cores)
        completed_commands += 1
        progress = (completed_commands / total_commands) * 100
        print(f"Progress: {progress:.2f}%")


if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Kmer Extraction Runner")
    parser.add_argument("-j", "--json", help="JSON file path")
    args = parser.parse_args()

    if args.json:
        main(args.json)
    else:
        parser.print_help()
