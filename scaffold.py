import numpy as np
import pandas as pd
import glob
import os
from tqdm import tqdm
import argparse
import amr_functions as amr


def main(file_path, output):

    #Confirms the creation of the output folder
    amr.create_folder(os.getcwd(), output)

    #Process fasta file into a scaffold for each isolate
    amr.scaffold_fasta_file(file_path, output)


if __name__ == '__main__':

    
    parser = argparse.ArgumentParser(description='Process the input FASTA file.')
    parser.add_argument("-f", "--filepath", default=os.path.join(os.getcwd(), 'raw_data\wgs\wgs-401.fasta'),
        help="Path to the input raw FASTA file (default: 'raw_data\wgs\wgs-401.fasta' in the current directory)")
    parser.add_argument("-o", "--output", default=os.path.join(os.getcwd(), 'scaffold'),
        help="Path to the scaffold output for each isolate (default: 'scaffold' in the current directory)")
    args = parser.parse_args()
    main(args.filepath, args.output)
