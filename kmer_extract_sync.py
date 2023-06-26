import numpy as np 
import pandas as pd 
import glob, os 
from tqdm import tqdm 
import argparse
import amr_functions as amr

def generate_kmers(seq, k):
    k_freq = {}
    for i in range(0, len(seq) - k + 1):
        kmer = seq[i:i + k]
        if kmer in k_freq:
            k_freq[kmer] += 1
        else:
            k_freq[kmer] = 1                
    return k_freq   

def process_files(k_mer, path, output, msg):
    """
    Generate k-mer frequencies for all files in a directory and store them in a csv file.

    Args:
        - k_mer (int): The value of k to use for generating k-mer frequencies.
        - path (str): The directory path containing the files to process.
        - output (str): The base file name to use for the output csv files.
        - msg (str): An additional message to add to the file name.

    Returns:
        None.
    """
    # Create the directories if they don't already exist
    output_for_file = os.path.join(output,msg)
    if not os.path.exists(output_for_file):
        os.makedirs(output_for_file)
        print(f"Folder '{msg}' created successfully in 'scaffold_genes'")
    else:
        print(f"Folder '{msg}' already exists in 'scaffold_genes'")



    tmp = os.getcwd()
    os.chdir(path)  # Set current working directory to given path
    file_names = []  # Create an empty list to store file names
    for k in range(k_mer, k_mer+1):  # Iterate over k-mer range
        full_df = pd.DataFrame()  # Create an empty pandas dataframe to store k-mer frequencies
        count = 0  # Initialize count variable to keep track of number of processed files
        for file in tqdm(glob.glob("*"), desc='Strain reading '):  # Iterate over files in current directory
            #print(file)  # Print the file name
            target_file = open(file)  # Open the file for reading
            file_names.append(file)  # Append the file name to the list
            df_row = pd.DataFrame(index=file_names)  # Create a new pandas dataframe with file names as index
            #next(target_file)  # Skip the first line of the file
            for line in target_file.readlines():  # Read the file line by line
                line = line.strip()  # Remove leading/trailing whitespace
                
                if line.startswith('>') or line.startswith('\n'):  # Skip lines starting with '>'
                    pass
                else:
        
                    #read_file = target_file.read()  # Read the remaining contents of the file
                    seq = "".join(line.split())  # Join the file contents and remove white spaces
                    kmer_freq = generate_kmers(seq, k)  # Generate k-mer frequencies for the sequence
                    kmer_freq = dict(sorted(kmer_freq.items()))  # Sort the dictionary of k-mer frequencies in alphabetical order
                    keys = full_df.keys()  # Get the column names of the dataframe
                    kmer_existed = []  # Create an empty list to store k-mers that already exist in the dataframe
                    for kmer in kmer_freq.keys():  # Iterate over k-mers in the current file
                        if kmer not in keys:  # If k-mer is not in the column names of the dataframe
                            kmer_existed.append(kmer)  # Append the k-mer to the list
                    df = pd.DataFrame(columns=kmer_existed)  # Create a new pandas dataframe with k-mer frequencies as columns
                    full_df = pd.concat([full_df,df], axis = 1)  # Concatenate the new dataframe with the existing dataframe along columns
            full_df = pd.concat([full_df, df_row], axis=0)  # Concatenate the new dataframe with the existing dataframe along rows
            
            for j in kmer_freq.keys():  # Iterate over k-mers in the current file
                full_df.at[file,j] = int(kmer_freq[j])  # Assign the k-mer frequency to the corresponding cell in the dataframe




            #print(full_df)
            target_file.close()  # Close the file
            count = count + 1  # Increment the count variable
            file_names.pop()  # Remove the last file name from the list
        os.chdir(tmp)
        #full_df.to_parquet(os.path.join(os.getcwd(), str(output) + str(msg) +str(k)+".parquet"))  # Write the dataframe to a parquet file with a specific file name format
        #print(full_df)
        full_df.to_csv(os.path.join(output_for_file, 'kmer'+str(k)+'.csv'), index=True, sep = ';')






parser = argparse.ArgumentParser()
parser.add_argument("-k", "--k_mer", type=int, choices=range(2, 12), default=5,
                    help="The k_mer argument (an integer from 2 to 11, default: 5)")
parser.add_argument("-p", "--path", default=os.path.join(os.getcwd(), 'scaffold_genes\\antibiotic_resistance'),
                    help="The path argument (default: 'scaffold' in the current directory)")
parser.add_argument("-o", "--output", default=os.path.join(os.getcwd(), 'lib\kmer\kmer_genes'),
                    help="The output path argument (default: 'lib_files\\kmer' in the current directory)")
parser.add_argument("-f", "--func",  default='__main__',
                        help="The function argument (default: __main__)")
parser.add_argument("-t", "--type",  default='amr',
                        help="The type argument (default: amr)")

args = parser.parse_args()
print("The value of the k_mer argument is:", args.k_mer)
print("The value of the path argument is:", args.path)
print("The value of the output path argument is:", args.output)
print("The value of the func argument is:", args.func)
print("The value of the type argument is:", args.type)
__name__ = args.func

if __name__ == "__main__":
    #python A_kmerExtract.py -f __main__ -k 5 -p scaffold_genes\antibiotic_resistance -o lib\kmer\kmer_genes -t amr
    #python A_kmerExtract.py -f __main__ -k 5 -p scaffold_genes\drug_target -o lib\kmer\kmer_genes -t dt
    #python A_kmerExtract.py -f __main__ -k 5 -p scaffold_genes\transporter -o lib\kmer\kmer_genes -t tpt
    #python A_kmerExtract.py -f __main__ -k 5 -p scaffold_genes\virulence_factor -o lib\kmer\kmer_genes -t vf
    process_files(args.k_mer, args.path,os.path.join(os.getcwd(),args.output), msg=args.type)

if __name__ == "__gene__":
    #python A_kmerExtract.py -f __gene__ -p lib\genes_amr -o scaffold_genes
    #python A_kmerExtract.py -f __main__ -p scaffold_genes -o lib\kmer\kmer_genes -k 3
    amr.parse_fasta(os.path.join(os.getcwd(),args.path),    
                    os.path.join(os.getcwd(),      
                    args.output))
