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
parser.add_argument("-k", "--k_mer", type=int, choices=range(2, 12), default=5,
                    help="The k_mer argument (an integer from 2 to 11, default: 5)")
parser.add_argument("-p", "--path", default=os.path.join(os.getcwd(), 'scaffold_genes\\antibiotic_resistance'),
                    help="The path argument (default: 'scaffold' in the current directory)")
parser.add_argument("-o", "--output", default=os.path.join(os.getcwd(), 'lib\kmer'),
                    help="The output path argument (default: 'lib_files\\kmer' in the current directory)")
parser.add_argument("-f", "--func",  default='main',
                        help="The function argument (default: __main__)")
parser.add_argument("-t", "--type",  default='amr',
                        help="The type argument (default: amr)")
parser.add_argument("-c", "--cores",  default='4',
                        help="The type number of cores used (default: 4)")

args = parser.parse_args()
print("The value of the k_mer argument is:", args.k_mer)
print("The value of the path argument is:", args.path)
print("The value of the output path argument is:", args.output)
print("The value of the func argument is:", args.func)
print("The value of the type argument is:", args.type)
func = args.func


if __name__ == '__main__':
    #python A_kmerExtract.py -f __main__ -k 5 -p scaffold_genes\antibiotic_resistance -o lib\kmer\kmer_genes -t amr
    #python A_kmerExtract.py -f __main__ -k 5 -p scaffold_genes\drug_target -o lib\kmer\kmer_genes -t dt
    #python A_kmerExtract.py -f __main__ -k 5 -p scaffold_genes\transporter -o lib\kmer\kmer_genes -t tpt
    #python A_kmerExtract.py -f __main__ -k 5 -p scaffold_genes\virulence_factor -o lib\kmer\kmer_genes -t vf
    if func == "kmer_mult":
        amr.generate_kmer_frequencies(k_mer=args.k_mer, 
                                      path=os.path.join(os.getcwd(),args.path),
                                      output=os.path.join(os.getcwd(),args.output), 
                                      msg=args.type, 
                                      threads=int(args.cores), 
                                      function_mult=amr.kmer_of_files_modular)
    if func == "gene-scaffold":
        #python A_kmerExtract.py -f __gene__ -p lib\genes_amr -o scaffold_genes
        #python A_kmerExtract.py -f __main__ -p scaffold_genes -o lib\kmer\kmer_genes -k 3
        amr.parse_fasta(fasta_folder=os.path.join(os.getcwd(),args.path),    
                        output_path=os.path.join(os.getcwd(),args.output))
    if func == "wgs-scaffold":
        amr.scaffold_fasta_file(file_path=os.path.join(os.getcwd(),args.path), 
                                output=os.path.join(os.getcwd(),args.output))
