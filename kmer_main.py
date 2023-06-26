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
                    help="The path argument (default: 'scaffold_genes' in the current directory)")
parser.add_argument("-o", "--output", default=os.path.join(os.getcwd(), 'lib\kmer'),
                    help="The output path argument (default: 'lib_files\\kmer' in the current directory)")
parser.add_argument("-f", "--func",  default='main-scaffold',
                        help="The function argument (default: main) (Op.: kmer_mult, kmer_sync, main-scaffold, gene-scaffold)")
parser.add_argument("-fn", "--foldername",  default='amr',
                        help="The Name of the folder to be created to store the output (for ngs-specialtygenes) (default: amr)")
parser.add_argument("-c", "--cores", type=int,  default='4',
                        help="The type number of cores used (default: 4)")

args = parser.parse_args()
print("The value of the k_mer argument is:", args.k_mer)
print("The value of the path argument is:", args.path)
print("The value of the output path argument is:", args.output)
print("The value of the func argument is:", args.func)
print("The value of the foldername argument is:", args.foldername)
func = args.func


if __name__ == '__main__':
    #python A_kmerExtract.py -f __main__ -k 5 -p scaffold_genes\antibiotic_resistance -o lib\kmer\kmer_genes -t amr
    #python A_kmerExtract.py -f __main__ -k 5 -p scaffold_genes\drug_target -o lib\kmer\kmer_genes -t dt
    #python A_kmerExtract.py -f __main__ -k 5 -p scaffold_genes\transporter -o lib\kmer\kmer_genes -t tpt
    #python A_kmerExtract.py -f __main__ -k 5 -p scaffold_genes\virulence_factor -o lib\kmer\kmer_genes -t vf

    #Generate k-mer frequencies for all files in a directory and store them in a CSV file.
    #Requires a modular function to multiprocess.
    if func == "kmer_mult": #Faster
        amr.generate_kmer_frequencies_mult(k_mer=args.k_mer, 
                                      path=os.path.join(os.getcwd(),args.path),
                                      output=os.path.join(os.getcwd(),args.output), 
                                      folder=args.type, 
                                      threads=int(args.cores), 
                                      function_mult=amr.kmer_of_files_modular)
        
    #Process a list of files and update a dataframe with k-mer frequencies. 
    #The modular aspect of this function is used to run in multiprocess exec.
    #When used without multiprocess it may take a lot longer to complete, specially if k_mer >9
    if func == "kmer_sync": #Slower
        # Get the current working directory and join it with the provided path
        directory = os.path.join(os.getcwd(), args.path)
        # Create an empty list to store file paths
        file_paths = []
        # Traverse the directory tree and append the paths of all files to file_paths list
        for root, _, files in os.walk(directory):
            file_paths.extend([os.path.join(root, file) for file in files])
        # Sort the file paths in ascending order
        file_paths.sort()
        # Generate sync DataFrame using the kmer_of_files_modular function
        sync_df = amr.kmer_of_files_modular(file_list=file_paths, dataframe=pd.DataFrame(), k=args.k_mer)
        # Save sync to a CSV file
        sync_df.to_csv(os.path.join(os.path.join(os.getcwd(),args.output), 'kmer' + str(args.k_mer) + '.csv'), index=True, sep=';')

    #Parse a FASTA file and split its sequences by strain name into separate output files.
    if func == "gene-scaffold":
        #python A_kmerExtract.py -f __gene__ -p lib\genes_amr -o scaffold_genes
        #python A_kmerExtract.py -f __main__ -p scaffold_genes -o lib\kmer\kmer_genes -k 3
        amr.parse_fasta(fasta_folder=os.path.join(os.getcwd(),args.path),    
                        output_path=os.path.join(os.getcwd(),args.output))
        
    #Processes a FASTA file into a scaffold of the isolate, extracts relevant information, and writes output to specified folder.
    if func == "wgs-scaffold":
        amr.scaffold_fasta_file(file_path=os.path.join(os.getcwd(),args.path), 
                                output=os.path.join(os.getcwd(),args.output))
