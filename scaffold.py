import numpy as np 
import pandas as pd 
import glob, os 
from tqdm import tqdm 
import argparse
import amr_functions as amr
import concurrent.futures
from multiprocessing import freeze_support
import multiprocessing



parser = argparse.ArgumentParser()
parser.add_argument("-p", "--path", default=os.path.join(os.getcwd(), 'raw_data\\ngs-specialtygenes'),
                    help="The path argument (default: 'scaffold_genes' in the current directory)")
parser.add_argument("-o", "--output", default=os.path.join(os.getcwd(), 'lib\\kmer'),
                    help="The output path argument (default: 'lib_files\\kmer' in the current directory)")
parser.add_argument("-f", "--func",  default='gene-scaffold',
                        help="The function argument (default: main) (Op.: kmer_mult, kmer_sync, main-scaffold, gene-scaffold)")


args = parser.parse_args()
func = args.func


if __name__ == '__main__':
    titles = {}
    titles = amr.read_ascii_art(os.path.join(os.getcwd(), 'titles.txt'))



    #Parse a FASTA file and split its sequences by strain name into separate output files.
    if func == "gene-scaffold":
        print(titles[func])
        #python A_kmerExtract.py -f __gene__ -p lib\genes_amr -o scaffold_genes
        #python A_kmerExtract.py -f __main__ -p scaffold_genes -o lib\kmer\kmer_genes -k 3
        amr.scaffold_gene(fasta_folder=os.path.join(os.getcwd(),args.path),    
                        output_path=os.path.join(os.getcwd(),args.output))
        
    #Processes a FASTA file into a scaffold of the isolate, extracts relevant information, and writes output to specified folder.
    if func == "wgs-scaffold":
        print(titles[func])
        #python kmer_main.py -f wgs-scaffold -p raw_data\ngs-wgs -o scaffold
        amr.scaffold_wgs(file_path=os.path.join(os.getcwd(),args.path), 
                                output=os.path.join(os.getcwd(),args.output))
