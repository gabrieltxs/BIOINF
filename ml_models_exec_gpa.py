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

import re


# set seed for reproducibility
os.environ["PL_GLOBAL_SEED"] = str(seed)
random.seed(seed)


parser = argparse.ArgumentParser()
#python E_MLmodelsExec.py
parser.add_argument("-k", "--k_mer", type=int, choices=range(2, 12), default=8,
                    help="The k_mer argument (an integer from 2 to 11, default: 5)")
'''parser.add_argument("-pk", "--pathk", default=os.path.join(os.getcwd(), 'lib\\kmer'),
                    help="The path kmerX argument (default: 'lib/kmer' in the current directory)")'''
parser.add_argument("-pk", "--pathk", default=os.path.join(os.getcwd(), 'lib\kmer'),
                    help="The path kmerX argument (default: 'lib/kmer' in the current directory)")
parser.add_argument("-pf", "--pathf", default=os.path.join(os.getcwd(), 'lib\\processed'),
                    help="The path features argument (default: 'lib/processed' in the current directory)")
parser.add_argument("-o", "--output", default=os.path.join(os.getcwd(), 'results'),
                    help="The output results path argument (default: 'results' in the current directory)")

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

# load list of selected samples
antibiotics = ['ceftazidime', 'ciprofloxacin','meropenem', 'tobramycin']

# List of folder names
folders = ['gexp']
#'amr', 'dt', 'tpt', 'vf', 'wgs'

############################# Variables initizalization ###################################

antibiotic_dfs = amr.load_antibiotic_dfs(input_file_path, antibiotics)

# call function to get list of models
models = amr.get_models()

header = 'index;ceftazidime;ciprofloxacin;meropenem;tobramycin;\n'

# define hyperparameter ranges
pbounds = {'n_estimators': (400, 800),
            'learning_rate': (0.01, 0.1),
            'max_depth': (3, 9),
            'num_leaves': (5, 50),
            'min_child_samples': (10, 100),
            'device':'gpu'}

############################# Escolha do K ###################################

# Iterates over the values of k to read the corresponding k-mer feature dataset
for folder in tqdm(folders):
    amr.create_folder(output_results_path, folder)
    # Reads the k-mer feature dataset file
    for model in models:
        # Prints the name of the current model being used
        print('\n\n'+type(model).__name__)
        
        # Opens the metadata file for the current model and delta values to store the results
        with open(os.path.join(output_results_path, folder, type(model).__name__ + str(args.k_mer) +'.txt'), 'a+') as f:
            try: 
                # Checks if the metadata file already exists
                with open(os.path.join(output_results_path, folder, type(model).__name__ + str(args.k_mer) +'.txt'), 'r') as f2:
                    # Reads the first line of the metadata file
                    primeiralinha = f2.readline()
                # Checks if the first line of the metadata file is the expected header, if not, writes it
                if primeiralinha != header:
                    f.write(header)
            except: 
                # Writes the header in case the metadata file does not exist
                f.write(header)
                
            f.write(f'\n{args.k_mer};')


            for antibiotic in tqdm(antibiotic_dfs):
                xis = pd.read_csv(os.path.join(input_kmer_path, folder ,'gexp.csv'), sep=',', header = 0)
                xis = xis.rename(columns = lambda x:re.sub('[^A-Za-z0-9_]+', '', x))
                # Sets the index of the DataFrame to the first column and drops it
                xis.set_index('Unnamed0', inplace=True, drop=True)
                xis.fillna(0, inplace=True)    


                filtered_xis = amr.filter_antibiotic_dfs(xis, antibiotic, antibiotic_dfs)
                

                antibiotic_dfs[antibiotic].fillna(0, inplace=True)
                yps = antibiotic_dfs[antibiotic][antibiotic]

                





                #amr.plot_classes(antibiotic_dfs[antibiotic],antibiotic,output_results_path)
                X_train, X_test, y_train, y_test = model_selection.train_test_split(filtered_xis, yps ,test_size=0.2, 
                                                                                    random_state=0, stratify=antibiotic_dfs[antibiotic][antibiotic])
                X_train = X_train.reset_index(drop=True)
                y_train = y_train.reset_index(drop=True)

                #sm = BorderlineSMOTE(random_state=42)
                #X_train, y_train = sm.fit_resample(X_train, y_train)
                #X_resampled, y_resampled = adasyn.fit_resample(X_train, y_train)
                full_df_08 = pd.concat([X_train, y_train], axis =1)
                amr.plot_classes(full_df_08, antibiotic, output_results_path)

                def lgbm_cv(n_estimators, learning_rate, max_depth, num_leaves, min_child_samples):
                    """
                    Run cross-validation on a training dataset using a LightGBM classifier with hyperparameters specified by the input.

                    Args:

                    - n_estimators (float): The number of boosting iterations.
                    - learning_rate: (float): The learning rate used in the boosting process.
                    - max_depth (float): The maximum depth of each decision tree in the ensemble.
                    - num_leaves (float): The maximum number of leaves in each decision tree.
                    - min_child_samples (float): The minimum number of samples required to be at a leaf node.
                    - X_train (pd.DataFrame): 
                    - y_train (pd.Series):

                    Returns:
                    - (float): The mean of macro F1 score for one fold of cross-validation using the specified hyperparameters.
                    """
                    # define hyperparameters as a dictionary
                    params = {'n_estimators': int(round(n_estimators)),
                            'learning_rate': learning_rate,
                            'max_depth': int(round(max_depth)),
                            'num_leaves': int(round(num_leaves)),
                            'min_child_samples': int(round(min_child_samples))}
                    
                    # return mean score from cross-validation
                    return np.mean(amr.run_cv(X_train, y_train, params))
                                
                pbounds = {'n_estimators': (400, 800),
                            'learning_rate': (0.01, 0.1),
                            'max_depth': (3, 9),
                            'num_leaves': (5, 50),
                            'min_child_samples': (10, 100)}

                # run Bayesian optimization
                optimizer = BayesianOptimization(f=lgbm_cv, pbounds=pbounds, random_state=42)
                optimizer.maximize(init_points=10, n_iter=2)

                # print best hyperparameters
                print(optimizer.max)
                f.write('%.3f;' % mean(optimizer.max['target']))

                space = amr.get_space_from_optimizer(optimizer)
                
                model = lgb.LGBMClassifier(**space, random_state=42)
                model.fit(X_train, y_train)

                y_pred = model.predict(X_test)
                score = f1_score(y_test, y_pred, average= 'weighted')  
                print('F1-Teste: '+str(score))
                f.write('%.3f;\n' % mean(score))
                f.write(classification_report(y_test, y_pred))

            f.write('\n\n')

