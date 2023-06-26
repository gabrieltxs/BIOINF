# Import necessary libraries for data processing
import glob, os
import random
import numpy as np
import pandas as pd
from tqdm import tqdm 
import seaborn as sns


from sklearn.preprocessing import MinMaxScaler

# Import libraries for machine learning
from sklearn.metrics import f1_score
from sklearn.utils import resample
from sklearn.model_selection import StratifiedKFold

# Import libraries for Bayesian optimization
import warnings
warnings.filterwarnings('ignore')

# Import classifier models
from lightgbm import LGBMClassifier
import lightgbm as lgb

# Import visualization libraries
import matplotlib.pyplot as plt

# Set seed for reproducibility
from random import seed
seed(42)

# set seed for reproducibility
os.environ["PL_GLOBAL_SEED"] = str(seed)
random.seed(seed)


def load_antibiotic_dfs(path, antibiotics):
    """
    Load antibiotic data from CSV files into a dictionary of Pandas DataFrames.

    Args:
    - path (str): the path to the directory containing the CSV files.
    - antibiotics (list of str): the names of the antibiotics to load.

    Returns:
    - antibiotic_dfs (dict of str:Pandas DataFrame): a dictionary of Pandas DataFrames containing the loaded data. The keys are the antibiotic names and the values are the corresponding DataFrames.
    """
    antibiotic_dfs = {} # Initialize an empty dictionary to store the loaded DataFrames.
    for antibiotic in antibiotics: # Iterate over the list of antibiotics to load.
        filename = os.path.join(path, antibiotic+'_AMR.csv') # Construct the filename of the CSV file.
        df = pd.read_csv(filename, index_col=False, sep=';', header=0) # Load the CSV file into a Pandas DataFrame.
        antibiotic_dfs[antibiotic] = df # Add the DataFrame to the dictionary using the antibiotic name as the key.
        print(antibiotic, df.shape) # Print the name of the antibiotic and the shape of the loaded DataFrame.
    return antibiotic_dfs # Return the dictionary of loaded DataFrames.

def balance_classes(df: pd.DataFrame, col: str, target_minority: float, target_majority: float) -> pd.DataFrame:
    """
    Balance the classes of a binary classification problem in a Pandas DataFrame by either upsampling the minority class or downsampling the majority class.

    Args:
    - df (Pandas DataFrame): the DataFrame containing the data to balance.
    - col (str): the name of the column containing the binary target variable.
    - target_minority (float): the target proportion of the minority class after balancing.
    - target_majority (float): the target proportion of the majority class after balancing.

    Returns:
    - df_balanced (Pandas DataFrame): the balanced DataFrame.
    """
    n_samples = len(df) # Calculate the total number of samples in the DataFrame.
    n_majority = int(n_samples * target_minority) # Calculate the target number of samples for the majority class.
    n_minority = int(n_samples * target_majority) # Calculate the target number of samples for the minority class.

    if (df[col] == 1).sum() > (df[col] == 0).sum():
        # If the majority class is labeled as 1 in the DataFrame, select the majority and minority classes accordingly.
        df_majority = df[df[col] == 1]
        df_minority = df[df[col] == 0]
    else:
        # Otherwise, select the majority and minority classes the other way around.
        df_majority = df[df[col] == 0]
        df_minority = df[df[col] == 1]

    # Upsample the minority class using the resample function from scikit-learn.
    df_minority_upsampled = resample(df_minority,
                                     replace=False,
                                     n_samples=n_minority,
                                     random_state=123)

    # Downsample the majority class using the resample function from scikit-learn.
    df_majority_downsampled = resample(df_majority,
                                       replace=True,
                                       n_samples=n_majority,
                                       random_state=123)

    if (df[col] == 1).sum() > (df[col] == 0).sum():
        # If the majority class is labeled as 1 in the DataFrame, concatenate the upsampled and downsampled DataFrames accordingly.
        df_balanced = pd.concat([df_majority_downsampled, df_minority_upsampled])
    else:
        # Otherwise, concatenate the upsampled and downsampled DataFrames the other way around.
        df_balanced = pd.concat([df_minority_upsampled, df_majority_downsampled])

    return df_balanced # Return the balanced DataFrame.

def re_numero(va: int) -> str:
    """
    Convert an integer variable representing the susceptibility of a bacteria strain to an antibiotic into a string variable with the corresponding label.

    Args:
    - va (int): the integer variable representing the susceptibility of the bacteria strain to the antibiotic. Should be either 1 (susceptible) or 0 (resistant).

    Returns:
    - label (str): the corresponding label for the susceptibility of the bacteria strain to the antibiotic. Either 'Susceptible' or 'Resistant'.
    """
    if va == 1:
        label = 'Susceptible' # If the input variable is 1, set the label to 'Susceptible'.
    else:
        label = 'Resistant' # Otherwise, set the label to 'Resistant'.
    return label # Return the label.

def plot_classes(df: pd.DataFrame, col: str, output_results_path):
    """
    Plot the distribution of classes in a pandas DataFrame for a specific column.

    Args:
    - df (pd.DataFrame): the pandas DataFrame to plot.
    - col (str): the column of the pandas DataFrame to plot.
    - output_results_path (str): output path
    - antibiotic (str): specific antibiotic

    Returns:
    - None
    """

    # EDA for the data
    counts = df[col].value_counts() # Count the frequency of each class in the selected column.
    total_counts = counts.sum() # Get the total number of observations in the column.

    # Set custom color scheme
    colors = ['#8c510a', '#01665e']

    # Set font sizes
    SMALL_SIZE = 14
    MEDIUM_SIZE = 16
    BIGGER_SIZE = 18

    # Set font family
    plt.rcParams["font.family"] = "Arial"

    # Create a new plot and set its size
    fig, ax = plt.subplots(figsize=(8,6))

    # Create a bar plot of the class counts using the custom color scheme
    ax.bar([counts.index[0],counts.index[1]], counts.values, color=colors)

    # Add annotations to the bars to display the exact counts and percentages
    for i, count in enumerate(counts.values):
        ax.annotate(f"{count}\n({count/total_counts*100:.1f}%)", xy=(counts.index[i], count), ha='center', va='bottom', fontsize=SMALL_SIZE)

    # Set the title and axis labels for the plot
    ax.set_title(f'{col} Resistance (n={total_counts})', fontsize=BIGGER_SIZE, pad=20)
    ax.set_xlabel('Resistance Category', fontsize=MEDIUM_SIZE)
    ax.set_ylabel('Count', fontsize=MEDIUM_SIZE)
    ax.set_xticks([counts.index[0], counts.index[1]])

    # Add captions to the labels
    ax.text(counts.index[0], -60, re_numero(counts.index[0]), fontsize=SMALL_SIZE, ha='center')
    ax.text(counts.index[1], -60, re_numero(counts.index[1]), fontsize=SMALL_SIZE, ha='center')

    # Set tick parameters
    ax.tick_params(axis='both', which='major', labelsize=SMALL_SIZE, length=8, width=2)

    # Set spine parameters
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Save the plot to the output file path
    plt.savefig(os.path.join(output_results_path, f'{col}_Plot_AMR.png'), dpi=300, bbox_inches='tight')
    plt.close()

def get_models():
	models = list()
	models.append(LGBMClassifier(max_depth=10, n_estimators=500))
	return models
    
def run_cv(X_train, y_train, params):
    """
    Run cross-validation on a training dataset using a LightGBM classifier.

    Args:
    - X_train (pandas.DataFrame): A pandas dataframe containing the training data features.
    - y_train (pandas.Series): A pandas series containing the training data target variable.
    - params (dict): A dictionary containing hyperparameters to be passed to the LightGBM classifier.

    Returns:
    - (list): A list of scores, where each score is the macro F1 score for one fold of cross-validation.
    """
   
    # define cross-validation splitter
    skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)

    # run cross-validation
    scores = []
    # use the StratifiedKFold object to split the data into training and validation sets for each fold
    for train_idx, val_idx in skf.split(X_train, y_train, groups=None):
        # use the indices to select the relevant subsets of the data
        train_X, val_X = X_train.loc[train_idx], X_train.loc[val_idx]
        train_y, val_y = y_train.loc[train_idx], y_train.loc[val_idx]

        # initialize a LightGBM classifier with the provided hyperparameters, and fit it to the training data
        model = lgb.LGBMClassifier(**params, random_state=42)
        model.fit(train_X, train_y, eval_set=[(val_X, val_y)], early_stopping_rounds=10, verbose=0)

        # make predictions on the validation set using the trained model, and compute the macro F1 score
        y_pred = model.predict(val_X)
        score = f1_score(val_y, y_pred, average='macro')
        scores.append(score)

    # return the list of scores
    return scores

def get_space_from_optimizer(optimizer):
    """
    Extract the hyperparameter space from a hyperparameter optimizer object.

    Parameters:
    - optimizer (skopt.optimizer.Optimizer): A skopt optimizer object containing the results of a Bayesian optimization.

    Returns:
    - (dict): A dictionary of hyperparameters that can be passed to a model training function.
    """
    # initialize an empty dictionary to hold the hyperparameters
    space = {}

    # iterate over the hyperparameters in the optimizer object's best set of parameters
    for item in optimizer.max['params']:
        # if the hyperparameter is not the learning rate, cast it to an integer and add it to the dictionary
        if item != 'learning_rate':
            space[item] = int(optimizer.max['params'][item])
        # otherwise, add it to the dictionary as is
        else:
            space[item] = optimizer.max['params'][item]
    
    # return the dictionary of hyperparameters
    return space

def filter_antibiotic_dfs(xis, antibiotic, antibiotic_dfs):
    """
    Filters a DataFrame xis based on a column in the antibiotic_dfs dictionary.

    Parameters:
    - xis (pd.DataFrame): DataFrame to filter.
    - antibiotic (str): Name of the antibiotic to filter on.
    - antibiotic_dfs (dict): Dictionary with DataFrames to filter on.

    Returns:
    - pd.DataFrame: Filtered DataFrame.
    """
    # Create a list with the genome names for the specified antibiotic
    genome_names = antibiotic_dfs[antibiotic]['Genome Name'].tolist()

    # Add ".fna" to the end of each item in the genome_names list
    genome_names = [name + '.fasta' for name in genome_names]

    # Filter the DataFrame to only include rows where the 'Genome Name' column is in the genome_names list
    xis = xis.loc[xis.index.isin(genome_names)]

    return xis

def create_folder(directory, folder):
    """
    Creates a  folder in the given directory if it does not already exist.

    Parameters:
        - directory (str): The directory path to check and create 'tmp' folder.
        - folder (str): The name of the folder

    Returns:
        None
    """
    path_tmp = os.path.join(directory, folder)
    if not os.path.exists(path_tmp):
        os.makedirs(path_tmp)
        print(f"Folder '{folder}' created in directory:", directory)
    else:
        print(f"Folder '{folder}' already exists in directory:", directory)

def parse_fasta(fasta_folder: str, output_path: str) -> None:
    """
    Parse a FASTA file and split its sequences by strain name into separate output files.

    Parameters:
    - fasta_folder (str): the path to the input FASTA folder
    - output_path (str): the path to the output directory

    Returns:
    - None
    """
    os.chdir(fasta_folder)

    for file in tqdm(glob.glob("*"), desc="Processing files", unit="file"):
        with open(file) as f:
            # Create the directories if they don't already exist
            file_name = file.split('.')[0]
            output_for_file = os.path.join(output_path, file_name)
            create_folder(output_path, file_name)

            current_strain = None
            current_sequence = ""
            total_lines = sum(1 for _ in f)  # Count total lines in the file
            f.seek(0)  # Reset file pointer to the beginning

            for line_idx, line in enumerate(tqdm(f, total=total_lines, desc="Reading lines", unit="line", leave=False)):
                line = line.strip()
                if line.startswith(">"):
                    # Get the strain name from the header line
                    current_strain = line.split('[Pseudomonas')[1].split('|')[0].strip()
                    current_strain = 'Pseudomonas ' + current_strain
                    # Open the output file for the current strain
                    output_file = os.path.join(output_for_file, f"{current_strain}.fasta")
                    if os.path.exists(output_file):
                        # If the file already exists, append to it
                        with open(output_file, "a") as out:
                            out.write(current_sequence)
                    else:
                        # If the file doesn't exist yet, create it and write to it
                        with open(output_file, "w") as out:
                            out.write(current_sequence)

                    # Reset the current sequence
                    current_sequence = '\n' + line + '\n'
                else:
                    # Append to the current sequence
                    current_sequence += line

                # Update progress every 10% of lines processed
                if line_idx % (total_lines // 10) == 0:
                    progress_percentage = (line_idx / total_lines) * 100
                    tqdm.write(f"Progress: {progress_percentage:.2f}%")

            # Write the last sequence to the current output file
            output_file = os.path.join(output_for_file, f"{current_strain}.fasta")
            if os.path.exists(output_file):
                with open(output_file, "a") as out:
                    out.write(current_sequence)
            else:
                with open(output_file, "w") as out:
                    out.write(current_sequence)

def add_df2_to_df1(df1, df2, output_path, k, multiplier = 1):
    """
    Adds the corresponding values from df2, multiplied by the specified multiplier, 
    to each cell in df1 at row 'xyz' and every column in df1 that also exists in df2.
    
    Args:
    df1: pandas dataframe, the dataframe to update
    df2: pandas dataframe, the dataframe to use for adding values to df1
    output_path: str, the path to save the updated df1 dataframe
    k: int, value of the size of k
    multiplier: int or float, the value to multiply each column in df2 by before adding it to df1
    
    Returns:
    None
    """
    # iterate over all columns in df1
    for row in tqdm(df2.index.tolist(), desc='Strains proccessed'):
        for col in df2.columns:

            # check if the specified column exists in both dataframes
            if col in df1.columns:
                # add every cell in df2 at row 'xyz', column 'col' to the equivalent cell in df1
                df1.loc[row, col] = df1.loc[row, col] * df2.loc[row, col] * multiplier
            #else:
                #print(f"The '{col}' column does not exist in df2.")
        
    # save the updated df1 dataframe to the specified output path
    df1.to_csv(os.path.join(os.getcwd(),output_path+'\\kmer'+str(k)+'.csv'), index=True, sep=';')

def count_files(directory):
    """
    Counts the number of files in a directory.

    Args:
        directory (str): Path to the directory.

    Returns:
        int: Number of files in the directory.
    """
    file_count = 0

    for _, _, files in os.walk(directory):
        file_count += len(files)

    return file_count

def divide_files_into_groups(directory, num_groups):
    """
    Divides the files in a directory into agglomerated groups.

    Args:
        directory (str): Path to the directory.
        num_groups (int): Number of groups to divide the files into.

    Returns:
        list: A list of lists containing the file paths in each group.
    """
    file_paths = []

    for root, _, files in os.walk(directory):
        file_paths.extend([os.path.join(root, file) for file in files])

    file_paths.sort()

    # Calculate the number of files per group
    files_per_group = len(file_paths) // num_groups

    # Divide the files into groups
    groups = []
    start_index = 0

    for i in range(num_groups - 1):
        end_index = start_index + files_per_group
        groups.append(file_paths[start_index:end_index])
        start_index = end_index

    # Add the remaining files to the last group
    groups.append(file_paths[start_index:])

    # Create empty dataframes for each group
    gp_dataframes = [pd.DataFrame() for _ in range(num_groups)]

    return groups, gp_dataframes

def generate_kmers(seq, k):
    """
    Generate k-mer frequencies for a given sequence.
    
    Args:
        seq (str): Input sequence.
        k (int): Length of the k-mer.
        
    Returns:
        dict: Dictionary containing k-mers as keys and their frequencies as values.
    """
    k_freq = {}  # Initialize an empty dictionary to store k-mer frequencies
    
    for i in range(0, len(seq) - k + 1):  # Iterate through the indices of the sequence to generate k-mers
        kmer = seq[i:i + k]  # Extract the current k-mer from the sequence

        if kmer in k_freq:  # If the k-mer is already present in the dictionary
            k_freq[kmer] += 1  # Increment its frequency by 1
        else:
            k_freq[kmer] = 1  # Add the k-mer to the dictionary with a frequency of 1

    return k_freq  # Return the dictionary of k-mer frequencies

def merge_dicts(dict1, dict2):
    """
    Merge two dictionaries and return the merged dictionary.
    
    Args:
        dict1 (dict): The first dictionary.
        dict2 (dict): The second dictionary.
        
    Returns:
        dict: The merged dictionary containing keys and values from both input dictionaries.
    """
    
    merged_dict = dict(dict1)  # Create a copy of the first dictionary
    
    for key, value in dict2.items():
        if key in merged_dict:
            merged_dict[key] += int(value)  # Add up the values if the key already exists
        else:
            merged_dict[key] = int(value)  # Insert new value if the key doesn't exist
    
    return merged_dict

def process_files_multithread(file_list, dataframe, k):
    """
    Process a list of files and update a dataframe with k-mer frequencies.
    
    Args:
        file_list (list): List of file paths to be processed.
        dataframe (pandas.DataFrame): Existing dataframe to be updated.
        k (int): k value of the chosen k-mer
        
    Returns:
        pandas.DataFrame: Updated dataframe with k-mer frequencies.
    """
    file_names = []  # Initialize an empty list to store file names
    full_df = dataframe  # Assign the input dataframe to a new variable
    
    for file in tqdm(file_list, desc='Strain reading '):  # Iterate through the files in the file list
        kmer_freq = {}
        target_file = open(file)  # Open the current file
        file_name = file.split('\\')[len(file.split('\\'))-1]
        file_names.append(file_name)  # Add the current file name to the list
        
        df_row = pd.DataFrame(index=file_names)  # Create a new dataframe row with the current file name as index

        for line in target_file.readlines():  # Iterate through each line in the current file
            line = line.strip()  # Remove leading and trailing whitespaces from the line

            if line.startswith('>') or line.startswith('\n'):  # If the line starts with '>' or is empty
                pass  # Skip the line and continue to the next iteration
            else:
                seq = "".join(line.split())  # Remove whitespaces within the line and concatenate the characters
                kmer_freq_new = generate_kmers(seq, k)  # Generate k-mer frequencies for the sequence using a defined function
                kmer_freq = merge_dicts(kmer_freq_new, kmer_freq)
                #print(kmer_freq[i] for i in range(0,9))
                kmer_freq = dict(sorted(kmer_freq.items()))  # Sort the k-mer frequencies dictionary by keys
                keys = full_df.keys()  # Get the column keys of the full dataframe
                kmer_existed = []  # Initialize an empty list to store existing k-mers

                for kmer in kmer_freq.keys():  # Iterate through the keys of the k-mer frequencies dictionary
                    if kmer not in keys:  # If the k-mer key is not present in the dataframe columns
                        kmer_existed.append(kmer)  # Add the k-mer to the list of existing k-mers

                df = pd.DataFrame(columns=kmer_existed)  # Create a new dataframe with columns for existing k-mers
                full_df = pd.concat([full_df, df], axis=1)  # Concatenate the new dataframe with the full dataframe

        full_df = pd.concat([full_df, df_row], axis=0)  # Concatenate the row dataframe with the full dataframe vertically

        for j in kmer_freq.keys():  # Iterate through the keys of the k-mer frequencies dictionary
            full_df.at[file_name, j] = int(kmer_freq[j])  # Update the corresponding cell in the full dataframe with the k-mer frequency

        target_file.close()  # Close the current file
        file_names.pop()  # Remove the current file name from the list

    return full_df  # Return the updated full dataframe

def normalize_dataframe(df):
    """
    Normalize the data in a dataframe using Min-Max scaling.

    Args:
        df (pandas.DataFrame): The input dataframe containing the data to be normalized.

    Returns:
        pandas.DataFrame: The normalized dataframe with scaled values.
    """
    # Create a MinMaxScaler object
    scaler = MinMaxScaler()

    # Extract the numerical columns from the dataframe
    numerical_cols = df.select_dtypes(include=['float64', 'int64']).columns

    # Apply Min-Max scaling to the numerical columns
    df[numerical_cols] = scaler.fit_transform(df[numerical_cols])

    # Return the normalized dataframe
    return df

def compare_and_add(df1, df2):
    """
    Compare two DataFrames and add their corresponding cells.
    
    Args:
        df1 (pandas.DataFrame): The first DataFrame.
        df2 (pandas.DataFrame): The second DataFrame.
        
    Returns:
        pandas.DataFrame: A new DataFrame with the sum of corresponding cells from df1 and df2.
    """
    
    # Verify if columns are the same
    if not df1.columns.equals(df2.columns):
        print("Columns are not the same.")
        return None

    # Verify if rows are the same
    if not df1.index.equals(df2.index):
        print("Rows are not the same.")
        return None

    # Add up the cells and create a new DataFrame
    df3 = pd.DataFrame(index=df1.index, columns=df1.columns)
    for index in tqdm(df1.index):
        for column in df1.columns:

            teste = np.nansum([df1.loc[index, column] , df2.loc[index, column]])
            df3.at[index,column] = teste

    return df3

def process_specialty_genes_data(filename, list_genomes_path):
    """
    Process SpecialtyGenes-WGS data.

    Args:
    filename (str): The name of the file containing the SpecialtyGenes-WGS data.
    list_genomes_path (str): the folder to get the names of the genomes

    Returns:
        list: A list of DataFrames, each DataFrame representing the gene counts for a specific property.
    """

    #######################################--- Preprocessing Data ---############################################################
    create_folder(os.path.join(os.getcwd(), 'lib'), 'genes_eda_data')
    output_genes_eda = os.path.join(os.getcwd(), 'lib', 'genes_eda_data')

    create_folder(os.path.join(os.getcwd(), 'lib\kmer'), 'gexp')
    df_input_ml = os.path.join(os.getcwd(), 'lib\kmer', 'gexp')

    # Load SpecialtyGenes-WGS into a dataframe
    specialtygenes_df = pd.read_csv(filename)
    # Filter empty values in the 'Gene' column
    specialtygenes_df = specialtygenes_df.dropna(subset=['Gene'])
    unique_genes = list(specialtygenes_df['Gene'].unique())
    # Get the top 10 most frequent genes
    top_genes = specialtygenes_df['Gene'].value_counts().head(10).index.tolist()
    # List of columns to keep
    columns_to_keep = ['Property', 'Genome Name', 'Gene', 'Product']
    # Filter the DataFrame to keep only the specified columns
    specialtygenes_df = specialtygenes_df[columns_to_keep]


    #########################################--- General EDA ---#####################################################

    # EDA for the 'Gene' column
    eda_text = "EDA for the 'Gene' column:\n"
    eda_text += "Total number of genes: {}\n".format(len(specialtygenes_df['Gene']))
    eda_text += "Unique genes: {}\n\n".format(specialtygenes_df['Gene'].nunique())

    # Plot the distribution of genes
    plt.figure(figsize=(8, 6))
    gene_counts = specialtygenes_df['Gene'].value_counts().head(10)
    gene_counts.plot(kind='bar')
    plt.title("Top 10 Most Frequent Genes")
    plt.xlabel("Gene")
    plt.ylabel("Count")
    plt.tight_layout()
    
    # Rotate x-axis labels by 45 degrees
    plt.tick_params(axis='x', rotation=45)

    # Display values on top of bars
    for i, count in enumerate(gene_counts):
        plt.text(i, count, str(count), ha='center', va='bottom')


    # Save the plot as PNG
    plt.savefig(os.path.join(output_genes_eda, 'gene_distribution.png'))
    plt.close()

    #########################################--- Specific EDA ---####################################################

    # Get the unique properties
    properties = specialtygenes_df.groupby('Property')['Gene'].groups.keys()

    # Create subplots
    fig, axs = plt.subplots(len(properties), figsize=(10, 5 * len(properties)))

    # Iterate over properties
    for i, property in enumerate(properties):
        # Get the top 10 most frequent genes for the current property
        top_genes = specialtygenes_df[specialtygenes_df['Property'] == property]['Gene'].value_counts().head(10)
        eda_text += f"Most frequent gene for each group in the {property} column:\n"
        eda_text += str(top_genes) + "\n"

        # Plot the bar chart for the current property
        ax = axs[i]
        ax.bar(top_genes.index, top_genes.values)
        # Rotate x-axis labels by 45 degrees
        ax.tick_params(axis='x', rotation=45)
        ax.set_title(f"Top 10 Genes for Property '{property}'")
        ax.set_xlabel("Gene")
        ax.set_ylabel("Count")

        # Rotate x-axis labels by 45 degrees
        ax.tick_params(axis='x', rotation=45)

        # Display values on top of bars
        for j, count in enumerate(top_genes.values):
            ax.text(j, count, str(count), ha='center', va='bottom')


    # Adjust spacing between subplots
    plt.tight_layout()

    # Save the plot as PNG
    plt.savefig(os.path.join(output_genes_eda, 'most_frequent_gene.png'))
    plt.close()

    ####################################--- Text EDA ---###############################################################

    # Save the EDA text to a text file
    with open(os.path.join(output_genes_eda, 'eda_text.txt'), 'w') as file:
        file.write(eda_text)

    ######################################--- Generate Dataframes for ML ---##########################################



    # Get the list of files in the directory \scaffold_genes\\antibiotic_resistance'
    file_list = os.listdir(list_genomes_path)

    
    df_return = []

    for propertie in properties:
        # Create a DataFrame from the list of data
        tmp_df = specialtygenes_df[specialtygenes_df['Property'] == propertie]
        df = pd.DataFrame(index=file_list, columns=unique_genes)
        # Iterate over each row
        for index in tqdm(df.index):
            tmp_df2 = tmp_df[tmp_df['Genome Name'] == index.split('.')[0]]
            list_of_genome_genes = tmp_df2['Gene'].values.tolist()
            # Iterate over each column
            for column in df.columns:
                #print(column)
                if column in list_of_genome_genes:
                    count = list_of_genome_genes.count(column)
                    df.at[index, column] = count
                        

        # Save the 'propertie' dataframe to CSV
        #print(df)
        df_return.append(df)
        df.to_csv(os.path.join(df_input_ml, f"{propertie}.csv"))


    merged_kmer = compare_and_add(df_return[0], df_return[1])
    merged_kmer = compare_and_add(merged_kmer, df_return[2])

    merged_kmer.to_csv(os.path.join(df_input_ml, "gexp.csv"))

    # Return the resulting dataframe
    return df_return

def patch_dataframe(results, main_df, main_columns):
    """
    Update the main DataFrame by adding missing columns and assigning k-mer frequencies.

    Args:
        results (list): List of DataFrames containing k-mer frequency data.
        main_df (pd.DataFrame): The main DataFrame to be updated.
        main_columns (list): List of columns in the main DataFrame.

    Returns:
        pd.DataFrame: The updated main DataFrame.
    """
    for df in results:
        # Get the missing columns in main_df
        missing_columns = [col for col in df.columns if col not in main_columns]
        
        # Add the missing columns to main_df with None as the default value
        main_df = main_df.reindex(columns=main_df.columns.union(missing_columns), fill_value=None)
        
        for line in df.index.tolist():
            for j in df.keys():
                # Iterate over k-mers in the current file
                value = df.loc[line, j]
                if not np.isnan(value):
                    # Assign the k-mer frequency to the corresponding cell in the DataFrame
                    main_df.at[line, j] = int(value)
    
    return main_df

def scaffold_fasta_file(file_path, output):
    """
    Processes a FASTA file into a scaffold of the isolate, extracts relevant information, and writes output to specified folder.

    Args:
        file_path (str): Path to the input FASTA file.
        output (str): Path to the output folder.

    Output:
        None
    """
    # Create the output folder if it doesn't exist
    create_folder(os.getcwd(), output)

    # Read the input FASTA file
    with open(file_path, 'r') as fasta:
        for item in tqdm(fasta):
            if item[0] == '>':
                # Extract information from the header line
                n_contig = item.split('contig_')[1].split()[0]
                fna_name = item.split('[')[1].split('|')[0]

                if n_contig == '1':
                    # Create or append to the output file for the first contig
                    output_file = open(os.path.join(output, fna_name.strip() + '.fasta'), 'a')
                    output_file.close()
            else:
                # Write sequence data to the respective output file
                with open(os.path.join(output, fna_name.strip() + '.fasta'), 'w') as output_file:
                    output_file.write(item)