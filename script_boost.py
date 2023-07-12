import subprocess
import json
import os
import argparse

def boost(kmer, pathk, output, func):
    command = f"python boost.py -k {kmer} -pk {pathk} -o {output} -f {func}"
    print(f"Running command: {command}")
    subprocess.call(command, shell=True)
    print(f"Command completed: {command}")


def main(json_file):
    json_path = os.path.join(os.getcwd(), json_file)
    with open(json_path) as json_data:
        data = json.load(json_data)

    # Variables to store the values
    kmer = None
    pathk = None
    output = None

    # Lists to store the values
    func_list = []

    # Assign values to variables or lists based on the context
    for key, value in data.items():
        if isinstance(value['value'], list):
            if key == 'func':
                func_list = value['value']
        else:
            if key == 'kmer':
                kmer = value['value']
            elif key == 'output':
                output = value['value']
            elif key == 'pathk':
                pathk = value['value']

    total_commands = len(func_list)
    completed_commands = 0

    for i, func in enumerate(func_list):
        boost(kmer, pathk, output, func)
        completed_commands += 1
        progress = (completed_commands / total_commands) * 100
        print(f"Progress: {progress:.2f}%")


if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Boost Runner")
    parser.add_argument("-j", "--json", help="JSON file path", default='start_data_boost.json')
    args = parser.parse_args()

    if args.json:
        main(args.json)
    else:
        parser.print_help()
