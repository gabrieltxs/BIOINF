import os
from itertools import combinations
import pandas as pd
import amr_functions as amr
from tqdm import tqdm
import argparse
import time

def get_name_folder(folders):
    """
    Get the name of the folder by joining its parts with an underscore.

    Args:
        folders (list): List of folder names.

    Returns:
        str: Joined folder name.
    """
    return '_'.join(folders)

def add_dataframes(dataframes):
    """
    Combine multiple dataframes into a single dataframe.

    Args:
        dataframes (list): List of dataframes.

    Returns:
        pandas.DataFrame: Combined dataframe.
    """
    result = None
    for df in dataframes:
        if result is None:
            result = df
        else:
            # Start the timer
            start_time = time.time()
           # print("\nAdding dataframes...")

            result =df.mul(result)            

            # Stop the timer
            end_time = time.time()
            # Calculate the elapsed time
            elapsed_time = end_time - start_time
            # Print the Adding
            #print(f"Dataframes Adding time: {elapsed_time:.2} seconds")
    return result

def get_folder_combinations(directory):
    """
    Get all combinations of folders within a directory.

    Args:
        directory (str): Directory path.

    Returns:
        list: List of folder combinations.
    """
    folder_list = [folder for folder in os.listdir(directory) if os.path.isdir(os.path.join(directory, folder)) and not folder.startswith(".git")]

    combinations_list = []
    for r in range(2, len(folder_list) + 1):
        combinations_list.extend(combinations(folder_list, r))
    return combinations_list

def main(input_kmer_path, kmer, output):
    """
    Main function to process the directory and generate boosted dataframes.

    Args:
        directory (str): Directory path to the type of kmers
        input_kmer_path (str): Input kmer path.
        kmer (int): Kmer value.
        output (str): Output path.
    """
    amr.create_folder(output, 'boost')
    path_boost = os.path.join(output, 'boost')
    combinations_list = get_folder_combinations(input_kmer_path)
    folder_list = [folder for folder in os.listdir(input_kmer_path) if os.path.isdir(os.path.join(input_kmer_path, folder)) and not folder.startswith(".git")]
    list_df = {}

     # Start the timer
    start_time = time.time()
    print("Loading dataframes...")

    for folder in folder_list:
        df = pd.read_parquet(os.path.join(input_kmer_path, folder, f'kmer{kmer}.parquet'))
        list_df[folder] = df

    # Stop the timer
    end_time = time.time()
    # Calculate the elapsed time
    elapsed_time = end_time - start_time
    # Print the result
    print(f"Loading dataframes time: {elapsed_time:.2} seconds")

    for combination in tqdm(combinations_list, desc='Processing combinations'):
        print(get_name_folder(combination))
        amr.create_folder(path_boost, get_name_folder(combination))
        result_df = pd.DataFrame()
        list_combination_df =[]
        for item in combination:
            list_combination_df.append(list_df[item])
        
        boosted_df = add_dataframes(list_combination_df)
        boosted_df.to_parquet(os.path.join(path_boost, str(get_name_folder(combination)), f'kmer{kmer}.parquet'), index=True)

if __name__ == '__main__':
    # Command-line argument parsing

    parser = argparse.ArgumentParser()
    parser.add_argument("-k", "--kmer", type=int, choices=range(2, 12), default=8,
                        help="The k_mer argument (an integer from 2 to 11, default: 5)")
    parser.add_argument("-pk", "--pathk", default=os.path.join(os.getcwd(), 'lib\\kmer'),
                        help="The path kmerX argument (default: 'lib/kmer' in the current directory)")
    parser.add_argument("-o", "--output", default=os.path.join(os.getcwd(), 'lib'),
                        help="The output results path argument (default: 'results' in the current directory)")

    args = parser.parse_args()
    
    print("The value of the kmer argument is:", args.kmer)
    print("The value of the path kmerx argument is:", args.pathk)
    print("The output results path argument is:", args.output)
  


    # Call the main function with parsed arguments
    main(input_kmer_path=args.pathk,
         kmer=args.kmer,
         output=args.output)
