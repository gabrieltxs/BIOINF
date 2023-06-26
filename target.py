import pandas as pd
import os
import matplotlib.pyplot as plt
from sqlalchemy import true

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output", default=os.path.join(os.getcwd(), 'lib\\processed'),
                    help="The Output path argument (default: 'lib/processed' in the current directory)")
parser.add_argument("-p", "--path", default=os.path.join(os.getcwd(), 'lib'),
                    help="The path argument input(.tsv) (default: 'lib' in the current directory)")

args = parser.parse_args()

print("The value of the output argument is:", args.output)
print("The value of the path argument is:", args.path)


# Construct paths to input and output files relative to the current working directory
input_file_path = os.path.join(args.path, 'AMRphenotype.tsv')
output_file_path = args.output

# Read data from the input file and store it in a DataFrame
df = pd.read_table(input_file_path, header=0)

# Create a dictionary to store DataFrames for each antibiotic
antibiotic_dfs = {}

# Iterate over unique antibiotic names in the DataFrame
for antibiotic in df['Antibiotic'].unique():
    # Create an empty DataFrame with specific column names
    antibiotic_df = pd.DataFrame(columns=['Genome Name', antibiotic])
    
    # Iterate over rows in the DataFrame for the current antibiotic
    for item in df[df['Antibiotic'] == antibiotic].index.values:
        # Determine the value based on the 'Resistant Phenotype'
        if df.loc[item]['Resistant Phenotype'] == 'Resistant':
            value = 0
        if df.loc[item]['Resistant Phenotype'] == 'Susceptible':
            value = 1 
        if df.loc[item]['Resistant Phenotype'] == 'Intermediate':
            value = -1
        
        # Add the value to the appropriate cell in the antibiotic_df DataFrame
        antibiotic_df = antibiotic_df.append({'Genome Name': df.loc[item]['Genome Name'], antibiotic: value}, ignore_index=True)
    
    # Drop rows with value -1
    antibiotic_df = antibiotic_df[antibiotic_df[antibiotic] != -1]
    
    # Add the antibiotic_df DataFrame to the antibiotic_dfs dictionary with the antibiotic name as the key
    antibiotic_dfs[antibiotic] = antibiotic_df
    
    # Save antibiotic_df to a CSV file
    #antibiotic_df.reset_index(inplace=True)
    antibiotic_df.to_csv(os.path.join(output_file_path, f'{antibiotic}_AMR.csv'), index=False, sep = ';')

    # EDA for the data
    counts = antibiotic_df[antibiotic].value_counts()
    total_counts = counts.sum()

    # Set custom color scheme
    colors = ['#8c510a', '#01665e']

    # Set font sizes
    SMALL_SIZE = 14
    MEDIUM_SIZE = 16
    BIGGER_SIZE = 18

    # Set font family
    plt.rcParams["font.family"] = "Arial"

    fig, ax = plt.subplots(figsize=(8,6))
    ax.bar(['0','1'], counts.values, color=colors)

    # Add annotations to the bars
    for i, count in enumerate(counts.values):
        ax.annotate(f"{count}\n({count/total_counts*100:.1f}%)", xy=(i, count), ha='center', va='bottom', fontsize=SMALL_SIZE)

    ax.set_title(f'{antibiotic} Resistance (n={total_counts})', fontsize=BIGGER_SIZE, pad=20)
    ax.set_xlabel('Resistance Category', fontsize=MEDIUM_SIZE)
    ax.set_ylabel('Count', fontsize=MEDIUM_SIZE)
    ax.set_xticks([0, 1])
    #ax.set_xticklabels(['Resistant', 'Susceptible'], fontsize=SMALL_SIZE)

    # Add captions to the labels
    ax.text(0, -60, 'Susceptible', fontsize=SMALL_SIZE, ha='center')
    ax.text(1, -60, 'Resistant', fontsize=SMALL_SIZE, ha='center')

    # Set tick parameters
    ax.tick_params(axis='both', which='major', labelsize=SMALL_SIZE, length=8, width=2)

    # Set spine parameters
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Save the plot to the output file path
    plt.savefig(os.path.join(output_file_path, f'{antibiotic}_AMR.png'), dpi=300, bbox_inches='tight')
    plt.close()