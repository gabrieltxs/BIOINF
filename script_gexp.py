import subprocess
import json
import os
import argparse


def run_gexp(filename, list_genomes):
    """
    Runs the gexp.py script with the provided filename and list of genomes.
    """
    command = f"python gexp.py -f {filename} -l {list_genomes}"
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
        if key == 'filename':
            filename = value['value']
        elif key == 'list_genomes':
            list_genomes = value['value']

    run_gexp(filename, list_genomes)


if __name__ == '__main__':
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Run gexp.py script with provided JSON file')
    parser.add_argument('-j', '--json', help='Path to the JSON file', required=True)
    args = parser.parse_args()

    # Call the main function with the JSON file path
    main(args.json)
