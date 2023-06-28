import subprocess
import json
import os
import argparse


def run_scaffold(main, file_path, output_path):
    command = f"python scaffold.py -f {main} -p {file_path} -o {output_path}"
    print(f"Running command: {command}")
    subprocess.call(command, shell=True)
    print(f"Command completed: {command}")


def main(json_file):
    json_path = os.path.join(os.getcwd(), json_file)
    with open(json_path) as json_data:
        data = json.load(json_data)

    # Assign values to variables or dictionaries based on the context
    for key, value in data.items():
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


if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Scaffold Runner")
    parser.add_argument("-j", "--json", help="JSON file path")
    args = parser.parse_args()

    if args.json:
        main(os.path.join(os.getcwd(),args.json))
    else:
        parser.print_help()