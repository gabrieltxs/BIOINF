import os
import re
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

base_path = os.path.join(os.getcwd(), 'results')  # base folder path

folders = ['amr', 'dt', 'tpt', 'vf', 'wgs', 'gexp']
folders = [folder for folder in os.listdir('results') if os.path.isdir(os.path.join('results', folder)) and not folder.startswith(".git")]
pattern = r'(\d+\.\d+);(\d+\.\d+);'
k = 8
# Create a figure and subplots for each folder
fig, axs = plt.subplots(len(folders), 1, figsize=(8, 6 * len(folders)), sharex=True)


# Define the columns for the dataframes
columns = ['ceftazidime', 'ciprofloxacin', 'meropenem', 'tobramycin']

# Create empty dataframes
validation_df = pd.DataFrame(columns=columns)
test_df = pd.DataFrame(columns=columns)

# Traverse the folders
for i, folder in enumerate(folders):
    folder_path = os.path.join(base_path, folder)

    validation_scores = []
    test_scores = []

    # Traverse the files in the folder
    for file_name in os.listdir(folder_path):
        file_path = os.path.join(folder_path, file_name)
        if re.search(str(k) + str('.txt'), file_name):

            # Open the file and extract the values using regular expressions
            with open(file_path, 'r') as file:
                contents = file.read()
                matches = re.findall(pattern, contents)
                if matches:
                    for j, match in enumerate(matches):
                        if j < 4:  # Limit to the first four matches
                            values = [float(value) for value in match]
                            validation_scores.append(values[0])
                            test_scores.append(values[1])

    # Append scores to respective dataframes
    validation_df.loc[folder] = validation_scores
    test_df.loc[folder] = test_scores

    # Save dataframes as CSV files
    validation_df.to_csv('lib\genes_eda_data\\validation_scores.csv', sep=';', index_label='Folder')
    test_df.to_csv('lib\genes_eda_data\\test_scores.csv', sep=';', index_label='Folder')
