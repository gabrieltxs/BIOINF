import subprocess
import json
import os
import argparse


def run_gexp(output, path):
    """
    Runs the target.py script with the provided output and path.
    """
    command = f"python target.py -o {output} -p {path}"
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
        if key == 'output':
            output = value['value']
        elif key == 'path':
            path = value['value']

    run_gexp(output, path)


if __name__ == '__main__':
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Run target.py script with provided JSON file')
    parser.add_argument('-j', '--json', help='Path to the JSON file', required=True)
    args = parser.parse_args()

    # Call the main function with the JSON file path
    main(args.json)
