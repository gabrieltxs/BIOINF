import numpy as np 
import pandas as pd 
import glob, os 
from tqdm import tqdm 
import argparse
import amr_functions as amr
import concurrent.futures
from multiprocessing import freeze_support
import multiprocessing


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')



parser = argparse.ArgumentParser()
parser.add_argument("-k", "--k_mer", type=int, choices=range(2, 12), default=5,
                    help="The k_mer argument (an integer from 2 to 11, default: 5)")
parser.add_argument("-p", "--path", default=os.path.join(os.getcwd(), 'scaffold\\antibiotic_resistance'),
                    help="The path argument (default: 'scaffold_genes' in the current directory)")
parser.add_argument("-o", "--output", default=os.path.join(os.getcwd(), 'lib\\kmer'),
                    help="The output path argument (default: 'lib_files\\kmer' in the current directory)")
parser.add_argument("-f", "--func",  default='kmer_mult_gene',
                        help="The function argument (default: main) (Op.: kmer_mult, kmer_sync)")
parser.add_argument("-fn", "--foldername",  default='amr',
                        help="The Name of the folder to be created to store the output (for ngs-specialtygenes) (default: amr)")
parser.add_argument("-c", "--cores", type=int,  default='16',
                        help="The number of cores used (default: 4)")


args = parser.parse_args()
func = args.func


if __name__ == '__main__':
    titles = {}
    titles = amr.read_ascii_art(os.path.join(os.getcwd(), 'titles.txt'))


    #Generate k-mer frequencies for all files in a directory and store them in a CSV file.
    #Requires a modular function to multiprocess.
    if func == "kmer_mult_gene": #Faster
        print(titles[func])
        amr.generate_kmer_frequencies_mult(k_mer=args.k_mer, 
                                      path=os.path.join(os.getcwd(),args.path),
                                      output=os.path.join(os.getcwd(),args.output), 
                                      folder=args.foldername, 
                                      threads=int(args.cores), 
                                      function_mult=amr.kmer_of_files_modular_genes)
        
    #Generate k-mer frequencies for all files in a directory and store them in a CSV file.
    #Requires a modular function to multiprocess.
    if func == "kmer_mult_wgs": #Faster
        print(titles[func])
        amr.generate_kmer_frequencies_mult(k_mer=args.k_mer, 
                                      path=os.path.join(os.getcwd(),args.path),
                                      output=os.path.join(os.getcwd(),args.output), 
                                      folder=args.foldername, 
                                      threads=int(args.cores), 
                                      function_mult=amr.kmer_of_files_modular_wgs)
        
    #Process a list of files and update a dataframe with k-mer frequencies. 
    #The modular aspect of this function is used to run in multiprocess exec.
    #When used without multiprocess it may take a lot longer to complete, specially if k_mer >9
    
    if func == "kmer_sync_gene": #Slower
        print(titles[func])
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
        sync_df = amr.kmer_of_files_modular_genes(file_list=file_paths, dataframe=pd.DataFrame(), k=int(args.k_mer), wgs = args.wgs)
        # Save sync to a CSV file
        sync_df.to_csv(os.path.join(os.path.join(os.getcwd(),args.output), 'kmer' + str(args.k_mer) + '.csv'), index=True, sep=';')

    if func == "kmer_sync_wgs": #Slower
        print(titles[func])
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
        sync_df = amr.kmer_of_files_modular_wgs(file_list=file_paths, dataframe=pd.DataFrame(), k=int(args.k_mer), wgs = args.wgs)
        # Save sync to a CSV file
        sync_df.to_csv(os.path.join(os.path.join(os.getcwd(),args.output), 'kmer' + str(args.k_mer) + '.csv'), index=True, sep=';')
