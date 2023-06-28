import os
import argparse
import amr_functions as amr


def main(filename, input_path):
    titles = {}
    titles = amr.read_ascii_art(os.path.join(os.getcwd(), 'titles.txt'))
    print(titles['gexp'])
    # Process specialty genes data
    result_df = amr.process_specialty_genes_data(filename, input_path)



if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Process specialty genes data')

    # Add optional arguments with default values
    parser.add_argument("-f", "--filename", type=str, default=os.path.join(os.getcwd(), 'raw_data\\gexp-specialtygenes\\specialty_genes.csv'),
                    help="Path to specialty genes CSV file, default: lib\specialty_genes\specialty_genes.csv)")
    parser.add_argument("-l", "--list_genomes", type=str, default=os.path.join(os.getcwd(), 'scaffold\\antibiotic_resistance'),
                    help="Path to the list of the genomes names, default: scaffold\\antibiotic_resistance)")

    args = parser.parse_args()

    # Call the main function with the parsed arguments
    main(args.filename, args.list_genomes)
