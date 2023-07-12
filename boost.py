from operator import index
import os  # Importing the 'os' module for operating system related functions
from itertools import combinations  # Importing the 'combinations' function from the 'itertools' module
import pandas as pd  # Importing the 'pandas' library for data manipulation
import amr_functions as amr  # Importing the 'amr_functions' module
from tqdm import tqdm  # Importing the 'tqdm' library for progress bars
import argparse  # Importing the 'argparse' module for command-line argument parsing
import time  # Importing the 'time' module for timing operations
from sklearn.preprocessing import MinMaxScaler


def get_name_folder(folders):
    """
    Get the name of the folder by joining its parts with an underscore.

    Args:
        folders (list): List of folder names.

    Returns:
        str: Joined folder name.
    """
    return '_'.join(folders)  # Joining the folder names in the list with underscores

def normalize_df(df):
    min_value = df.min().min()
    max_value = df.max().max()
    normalized_df = (df - min_value) / (max_value - min_value)
    #normalized_df = normalized_df.round(3)
    return normalized_df


def add_dataframes(dataframes):
    """
    Combine multiple dataframes into a single dataframe.

    Args:
        dataframes (list): List of dataframes.

    Returns:
        pandas.DataFrame: Combined dataframe.
    """
    result = dataframes[0]  # Assign the first DataFrame as the initial result
    for df in dataframes[1:]:
        # Start the timer
        start_time = time.time()  # Recording the start time of the operation

        result = df.mul(result)  # Multiplying the current dataframe with the result dataframe

        # Stop the timer
        end_time = time.time()  # Recording the end time of the operation
        # Calculate the elapsed time
        elapsed_time = end_time - start_time  # Calculating the elapsed time
        # Print the Adding
        # print(f"Dataframes Adding time: {elapsed_time:.2} seconds")  # Printing the elapsed time for the operation

    return result  # Returning the combined dataframe


def get_folder_combinations(directory):
    """
    Get all combinations of folders within a directory.

    Args:
        directory (str): Directory path.

    Returns:
        list: List of folder combinations.
    """
    folder_list = [
        folder for folder in os.listdir(directory) if os.path.isdir(os.path.join(directory, folder)) and not folder.startswith(".git")
    ]  # Getting the list of folders in the input kmer path

    combinations_list = []
    for r in range(2, len(folder_list) + 1):
        combinations_list.extend(combinations(folder_list, r))  # Generating all combinations of folders
    return combinations_list  # Returning the list of folder combinations


def mult(input_kmer_path, kmer, output):
    """
    Main function to process the directory and generate boosted dataframes.

    Args:
        input_kmer_path (str): Input kmer path.
        kmer (int): Kmer value.
        output (str): Output path.
    """

    # Creating a 'boost' folder in the output directory
    amr.create_folder(output, 'boost')  
    # Creating the full path to the 'boost' folder
    path_boost = os.path.join(output, 'boost')  
    # Getting the list of folder combinations
    combinations_list = get_folder_combinations(input_kmer_path) 
    # Getting the list of folders in the input kmer path
    folder_list = [folder for folder in os.listdir(input_kmer_path) if os.path.isdir(os.path.join(input_kmer_path, folder)) and not folder.startswith(".git")]
    # Creating an empty dictionary to store the dataframes
    list_df = {}  

    for folder in tqdm(folder_list,desc='Loading dataframes'):
        # Reading the parquet file into a dataframe
        df = pd.read_parquet(os.path.join(input_kmer_path, folder, f'kmer{kmer}.parquet'))
        # Storing the dataframe in the dictionary with the folder name as the key
        list_df[folder] = df



    for combination in tqdm(combinations_list, desc='Processing combinations'):
        # Printing the name of the current combination
        print(get_name_folder(combination))
        # Creating a folder for the current combination in the 'boost' folder
        amr.create_folder(path_boost, get_name_folder(combination))
        # Creating an empty dataframe
        result_df = pd.DataFrame()
        # Creating an empty list to store the dataframes in the current combination
        list_combination_df = []
        for item in combination:
            # Appending the dataframes corresponding to the current combination to the list
            list_combination_df.append(list_df[item])

        # Combining the dataframes in the current combination into a single dataframe
        boosted_df = add_dataframes(list_combination_df)
        # Saving the boosted dataframe as a parquet file
        boosted_df.to_parquet(os.path.join(path_boost, str(get_name_folder(combination)), f'kmer{kmer}.parquet'), index=True)


def bond(kmer, output):
    """
    Main function to process the directory and generate boosted dataframes.

    Args:
        input_kmer_path (str): Input kmer path.
        kmer (int): Kmer value.
        output (str): Output path.
    """

    # Step 1: Load CSV file into 'gexp' DataFrame
    #csv_path = os.path.join(input_kmer_path,'gexp')
    gexp = pd.read_csv(os.path.join(output,'gexp','gexp.csv'), sep=';', header=0)
    gexp.set_index('Unnamed: 0', inplace=True, drop=True)    
    # Step 2: Iterate over list of directory names
    directory_list = os.listdir(os.path.join(output, 'boost'))

    # Create an empty dictionary to store the DataFrames
    dfs = {}

    amr.create_folder(output,"bond")
    bond_path = os.path.join(output, "bond")
    # Iterate over each directory
    for directory_name in tqdm(directory_list, desc='Loading top performers (boost)'):
        # Generate the file path for each Parquet file using the directory name
        parquet_path = os.path.join(output,'boost',directory_name, f'kmer{kmer}.parquet')  # Replace 'k' with the appropriate value
        
        # Step 3: Load Parquet file into DataFrame
        df = pd.read_parquet(parquet_path)
        
        # Store the DataFrame in the dictionary with the directory_name as the key
        dfs[directory_name] = df

    for df in tqdm(dfs.keys(), desc='Bonding'):
        # Create a MinMaxScaler object
        scaler = MinMaxScaler()

        # Fit the scaler on the gexp DataFrame and transform it
        #gexp_normalized = pd.DataFrame(scaler.fit_transform(gexp), columns=gexp.columns, index=gexp.index)
        gexp_normalized = normalize_df(gexp)
        # Fit the scaler on the dfs[df] DataFrame and transform it
        #dfs_normalized = pd.DataFrame(scaler.fit_transform(dfs[df]), columns=dfs[df].columns, index=dfs[df].index)
        dfs_normalized = normalize_df(dfs[df])                                                     
        # Step 4: Concatenate DataFrames horizontally
        concatenated_df = pd.concat([gexp_normalized, dfs_normalized], axis=1)


        # Generate the file path for the concatenated Parquet file
        amr.create_folder(bond_path, get_name_folder([df, "gexp"]))
        output_path = os.path.join(bond_path,get_name_folder([df, "gexp"]), f'kmer{kmer}.parquet')

        # Step 6: Save the concatenated DataFrame as a Parquet file
        concatenated_df.to_parquet(output_path)
        #concatenated_df.to_csv(output_path, index=True, sep=';')

        # Optional: Print the concatenated DataFrame
        print(concatenated_df)

    pass

        


if __name__ == '__main__':
    # Command-line argument parsing

    parser = argparse.ArgumentParser()
    parser.add_argument("-k", "--kmer", type=int, choices=range(2, 12), default=8,
                        help="The k_mer argument (an integer from 2 to 11, default: 5)")
    parser.add_argument("-pk", "--pathk", default=os.path.join(os.getcwd(), 'lib\\kmer'),
                        help="The path kmerX argument (default: 'lib/kmer' in the current directory)")
    parser.add_argument("-o", "--output", default=os.path.join(os.getcwd(), 'lib'),
                        help="The output results path argument (default: 'results' in the current directory)")
    parser.add_argument("-f", "--func", default='mult',
                        help="The type of function the boost.py shall perform boost (Mult, Bond) default: 'Mult', Bond")

    args = parser.parse_args()
    
    print("The value of the kmer argument is:", args.kmer)
    print("The value of the path kmerx argument is:", args.pathk)
    print("The output results path argument is:", args.output)
    print("The func argument is:", args.func)

    if args.func =='mult':
        titles = {}
        titles = amr.read_ascii_art(os.path.join(os.getcwd(), 'titles.txt'))
        print(titles['boosttimes'])

        # Call the mult function with parsed arguments
        mult(input_kmer_path=amr.get_absolute_path(args.pathk),
            kmer=args.kmer,
            output=amr.get_absolute_path(args.output))
        
    elif args.func =='bond':
        titles = {}
        titles = amr.read_ascii_art(os.path.join(os.getcwd(), 'titles.txt'))
        print(titles['bond'])
        # Call the mult function with parsed arguments
        bond(kmer=args.kmer,
            output=amr.get_absolute_path(args.output))