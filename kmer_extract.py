import numpy as np 
import pandas as pd 
import glob, os 
from tqdm import tqdm 
import argparse
import amr_functions as amr
import concurrent.futures
from multiprocessing import freeze_support
import multiprocessing


import multiprocessing

def main(k_mer, path, output, msg, threads):
    """
    Generate k-mer frequencies for all files in a directory and store them in a csv file.

    Args:
        - k_mer (int): The value of k to use for generating k-mer frequencies.
        - path (str): The directory path containing the files to process.
        - output (str): The base file name to use for the output csv files.
        - msg (str): An additional message to add to the file name.
        - threads (int): Number of processes for multiprocessing

    Returns:
        None.
    """
    # Create the directories if they don't already exist
    output_for_file = os.path.join(output,msg)
    amr.create_folder(output, msg)

    file_groups, df_groups = amr.divide_files_into_groups(path, threads)

    # Create a multiprocessing Pool with the desired number of processes
    with multiprocessing.Pool(processes=threads) as pool:
        # Map the tasks to the pool of processes
        results = pool.starmap(amr.process_files_multithread, zip(file_groups, df_groups, [k_mer]*threads))

    # Create the ultron_df dataframe
    ultron_df = pd.DataFrame()
    # Store the columns of ultron_df in a set for faster membership checking
    ultron_columns = set(ultron_df.columns)
    
    ultron_df = amr.patch_dataframe(results, ultron_df, ultron_columns)


    ultron_df.to_csv(os.path.join(output_for_file, 'kmer'+str(k_mer)+'.csv'), index=True, sep = ';')    


        





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
    if func == "main":
        main(args.k_mer, args.path,os.path.join(os.getcwd(),args.output), msg=args.type, threads=int(args.cores))
    if func == "gene":
        #python A_kmerExtract.py -f __gene__ -p lib\genes_amr -o scaffold_genes
        #python A_kmerExtract.py -f __main__ -p scaffold_genes -o lib\kmer\kmer_genes -k 3
        amr.parse_fasta(os.path.join(os.getcwd(),args.path),    
                        os.path.join(os.getcwd(),      
                        args.output))
