
import re
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
#import A_kmerExtract as k_extract

# Set seed for reproducibility
from random import seed
seed(42)

# Parse command-line arguments
import argparse


# set seed for reproducibility
os.environ["PL_GLOBAL_SEED"] = str(seed)
random.seed(seed)


parser = argparse.ArgumentParser()
parser.add_argument("-k", "--k_mer", type=int, choices=range(2, 12), default=8,
                    help="The k_mer argument (an integer from 2 to 11, default: 5)")
parser.add_argument("-pk", "--pathk", default=os.path.join(os.getcwd(), 'lib\\kmer'),
                    help="The path kmerX argument (default: 'lib/kmer' in the current directory)")
parser.add_argument("-o", "--output", default=os.path.join(os.getcwd(), 'results'),
                    help="The output results path argument (default: 'results' in the current directory)")
parser.add_argument("-pg", "--pathg", default=os.path.join(os.getcwd(), 'lib\\genes_amr'),
                    help="The genes_amr path argument (default: 'lib\\gene_amr' in the current directory)")


args = parser.parse_args()

print("The value of the k_mer argument is:", args.k_mer)
print("The value of the path kmerx argument is:", args.pathk)
print("The output results path argument is:", args.output)
print("The genes_amr path argument is:", args.pathg)
 
 

# Construct paths to input and output files relative to the current working directory
kmer_value = args.k_mer
input_kmer_path = args.pathk
output_results_path = args.output
input_genes_kmer = args.pathg

# load list of selected samples
antibiotic_dfs ={}
antibiotics = ['ceftazidime', 'ciprofloxacin','meropenem', 'tobramycin']


path_tmp=os.path.join(os.getcwd(), 'lib\\tmp')
amr.create_folder(os.path.join(os.getcwd(), 'lib'), 'tmp')

#for antibiotic in antibiotics:
    #antibiotic_dfs[antibiotic] = pd.read_csv(os.path.join(input_genes_kmer,antibiotic+'.csv'), index_col=False, sep=';', header = 0)

# Iterates over the values of k to read the corresponding k-mer feature dataset
for i in tqdm(range(args.k_mer,args.k_mer+1), desc='Boosting kmer'):
    kmer_read_amr ={}
    kmer =   pd.read_csv(os.path.join(input_kmer_path,'kmer' +str(i) +'.csv'), sep=';', header = 0)
    kmer.set_index('Unnamed: 0', inplace=True)
    print(kmer)
    kmer_gene = pd.read_csv(os.path.join(input_kmer_path, 'kmer_genes\\kmer_amr'+str(i) +'.csv'), sep=';', header =0)
    kmer_gene.set_index('Unnamed: 0', inplace=True)
    print(kmer_gene)
    amr.add_df2_to_df1(kmer, kmer_gene, output_results_path, i, 1)




