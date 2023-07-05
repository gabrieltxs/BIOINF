import subprocess
import json
import os
import argparse


def run_models(kmer, pathk, pathf, output, entry):
    """
    Runs the gexp.py script with the provided filename and list of genomes.
    """
    command = f"python models_exec.py -k {kmer} -pk {pathk} -pf {pathf} -o {output} -e {entry}"
    print(f"Running command: {command}")
    subprocess.call(command, shell=True)
    print(f"Command completed: {command}")


def main(json_file):
    """
    Main function that loads data from the JSON file and executes the script.
    """
    json_path = os.path.join(os.getcwd(), json_file)

    with open(json_path) as json_data:
        data = json.load(json_data)

    # Assign values to variables or dictionaries based on the context
    for key, value in data.items():
        if isinstance(value['value'], list):
            if key == 'entry':
                entry_list = value['value']
        else:
            if key == 'kmer':
                kmer = value['value']
            if key == 'pathk':
                pathk = value['value']
            if key == 'pathf':
                pathf = value['value']
            if key == 'output':
                output = value['value']

    for entry in entry_list:
        
        run_models(kmer, pathk, pathf, output, entry)


if __name__ == '__main__':
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Run models_exec.py script with provided JSON file')
    parser.add_argument('-j', '--json', help='Path to the JSON file', required=True)
    args = parser.parse_args()

    # Call the main function with the JSON file path
    main(args.json)
