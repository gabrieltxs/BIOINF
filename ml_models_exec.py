# Import necessary libraries for data processing
import os
import random
from tqdm import tqdm
import numpy as np
import pandas as pd
from numpy import mean
from imblearn.over_sampling import BorderlineSMOTE
# Import libraries for machine learning
from sklearn import model_selection
from sklearn.metrics import f1_score, classification_report
from sklearn.model_selection import StratifiedKFold
# Import libraries for Bayesian optimization
import warnings
warnings.filterwarnings('ignore')
from bayes_opt import BayesianOptimization
# Import classifier models
import lightgbm as lgb
# Import visualization libraries
import matplotlib.pyplot as plt
# Import amr_functions libraries
import amr_functions as amr
# Set seed for reproducibility
from random import seed
seed(42)
# Parse command-line arguments
import argparse


# set seed for reproducibility
os.environ["PL_GLOBAL_SEED"] = str(seed)
random.seed(seed)


parser = argparse.ArgumentParser()
#python E_MLmodelsExec.py
parser.add_argument("-k", "--k_mer", type=int, choices=range(2, 12), default=5,
                    help="The k_mer argument (an integer from 2 to 11, default: 5)")
'''parser.add_argument("-pk", "--pathk", default=os.path.join(os.getcwd(), 'lib\\kmer'),
                    help="The path kmerX argument (default: 'lib/kmer' in the current directory)")'''
parser.add_argument("-pk", "--pathk", default=os.path.join(os.getcwd(), 'lib'),
                    help="The path kmerX argument (default: 'lib/kmer' in the current directory)")
parser.add_argument("-pf", "--pathf", default=os.path.join(os.getcwd(), 'lib\\target'),
                    help="The path features argument (default: 'lib/target' in the current directory)")
parser.add_argument("-o", "--output", default=os.path.join(os.getcwd(), 'results'),
                    help="The output results path argument (default: 'results' in the current directory)")
parser.add_argument("-e", "--entry", default='gexp',
                    help="The argument that defines the type of entry of the ml model (default = 'kmer')")

#python E_MLmodelsExec.py -k 5 -pk lib\kmer -pf lib\\processed -o results
args = parser.parse_args()

print("The value of the k_mer argument is:", args.k_mer)
print("The value of the path kmerx argument is:", args.pathk)
print("The value of the path feature argument is:", args.pathf)
print("The output results path argument is:", args.output)


# Construct paths to input and output files relative to the current working directory
input_file_path = args.pathf
input_kmer_path = args.pathk
output_results_path = args.output
entry = args.entry
# load list of selected samples
antibiotics = ['ceftazidime', 'ciprofloxacin','meropenem', 'tobramycin']

# List of folder names
folders = amr.get_folders(input_kmer_path)
#, 'wgs'



############################# Variables initizalization ###################################

antibiotic_dfs = amr.load_antibiotic_dfs(input_file_path, antibiotics)

# call function to get list of models
models = amr.get_models()

header = 'index;ceftazidime;ciprofloxacin;meropenem;tobramycin;\n'

############################# Escolha do K ###################################

# Iterates over the values of k to read the corresponding k-mer feature dataset
for folder in tqdm(folders):
    amr.create_folder(output_results_path, folder)
    # Reads the k-mer feature dataset file
    for model in models:
        # Prints the name of the current model being used
        print('\n\n'+type(model).__name__)
        amr.process_model_results(output_results_path=output_results_path,
                                  folder=folder,
                                  model=model,
                                  kmer=args.k_mer,
                                  antibiotic_dfs=antibiotic_dfs,
                                  input_kmer_path=input_kmer_path,
                                  header=header,
                                  entry=entry)