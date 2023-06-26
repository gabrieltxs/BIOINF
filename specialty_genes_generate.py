import os
import argparse
import amr_functions as amr


def main(filename, input_path):

    # Process specialty genes data
    result_df = amr.process_specialty_genes_data(filename, input_path)

    # Print the resulting DataFrame
    print("result_df:")
    print(result_df)


if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Process specialty genes data')

    # Add optional arguments with default values
    parser.add_argument("-f", "--filename", type=str, default=os.path.join(os.getcwd(), 'lib\specialty_genes\specialty_genes.csv'),
                    help="Path to specialty genes CSV file, default: lib\specialty_genes\specialty_genes.csv)")
    parser.add_argument("-i", "--input_path", type=str, default=os.path.join(os.getcwd(), 'scaffold_genes\\antibiotic_resistance'),
                    help="Path to specialty genes CSV file, default: scaffold_genes\\antibiotic_resistance)")

    args = parser.parse_args()

    # Call the main function with the parsed arguments
    main(args.filename, args.input_path)
