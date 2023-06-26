import numpy as np
import pandas as pd
import glob
import os
from tqdm import tqdm
import argparse
import amr_functions as amr


def main(file_name):
    path = os.path.join('D:\\OneDrive\\Documentos\\OneDrive\\Documentos\\Trabalho de Conclusão de Curso\\DATASET\\', file_name)

    amr.create_folder(os.getcwd(), 'scaffold')

    if not os.path.exists('scaffold'):
        os.mkdir('scaffold')

    with open(path, 'r') as fasta, open('scaffold\\Pseudomonas aeruginosa strain CH4443.fna', 'w') as output_file:
        for item in tqdm(fasta):
            if item[0] == '>':
                n_contig = item.split('contig_')[1].split()[0]

                if n_contig == '1':
                    output_file.close()

                    fna_name = item.split('[')[1].split('|')[0]
                    output_file = open(os.path.join('scaffold', fna_name.strip() + '.fna'), 'w')
            else:
                output_file.write(item)



if __name__ == '__main__':

    

    parser = argparse.ArgumentParser(description='Process the input FASTA file.')
    parser.add_argument("-f", "--filename", default=os.path.join(os.getcwd(), 'results'),
        help="The output results path argument (default: 'results' in the current directory)")
    parser.add_argument('file_name', type=str, help='Name of the input FASTA file')
    args = parser.parse_args()
    main(args.file_name)
