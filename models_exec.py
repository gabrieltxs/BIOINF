import os
import random
from tqdm import tqdm
import numpy as np
import pandas as pd
from numpy import mean
from imblearn.over_sampling import BorderlineSMOTE
from sklearn import model_selection
from sklearn.metrics import f1_score, classification_report
from sklearn.model_selection import StratifiedKFold
import warnings
warnings.filterwarnings('ignore')
from bayes_opt import BayesianOptimization
import lightgbm as lgb
import matplotlib.pyplot as plt
import amr_functions as amr
from random import seed
seed(42)
import argparse

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-k", "--k_mer", type=int, choices=range(2, 12), default=8,
                        help="The k_mer argument (an integer from 2 to 11, default: 5)")
    parser.add_argument("-pk", "--pathk", default=os.path.join(os.getcwd(), 'lib'),
                        help="The path kmerX argument (default: 'lib/kmer' in the current directory)")
    parser.add_argument("-pf", "--pathf", default=os.path.join(os.getcwd(), 'lib\\target'),
                        help="The path features argument (default: 'lib/target' in the current directory)")
    parser.add_argument("-o", "--output", default=os.path.join(os.getcwd(), 'results'),
                        help="The output results path argument (default: 'results' in the current directory)")
    parser.add_argument("-e", "--entry", default='bond',
                        help="The argument that defines the type of entry of the ml model (default = 'kmer', 'boost')")
    
    args = parser.parse_args()
    
    print("The value of the k_mer argument is:", args.k_mer)
    print("The value of the path kmerx argument is:", args.pathk)
    print("The value of the path feature argument is:", args.pathf)
    print("The output results path argument is:", args.output)
    
    input_file_path = args.pathf
    input_kmer_path = args.pathk
    output_results_path = args.output
    entry = args.entry
    
    os.environ["PL_GLOBAL_SEED"] = str(seed)
    random.seed(seed)
    
    # Load list of selected samples
    antibiotics = amr.extract_antibiotic_names(input_file_path)
    folders = amr.get_folders(os.path.join(input_kmer_path,entry))
    antibiotic_dfs = amr.load_antibiotic_dfs(input_file_path, antibiotics)
    models = amr.get_models()



    if entry == 'kmer':
        for folder in tqdm(folders):
            amr.create_folder(output_results_path, folder)
            titles = {}
            titles = amr.read_ascii_art(os.path.join(os.getcwd(), 'titles.txt'))
            print(titles['models_exec_kmer'])
            for model in models:
                for k in range(8, 9):
                        print(f"Kmer Analysis: {k}")
                        print('\n\n'+type(model).__name__)
                        amr.process_model_results(output_results_path=output_results_path,
                                                folder=folder,
                                                model=model,
                                                kmer=k,
                                                antibiotic_dfs=antibiotic_dfs,
                                                input_kmer_path=input_kmer_path,
                                                entry=entry)
    if entry == 'boost':
        for folder in tqdm(folders):
            amr.create_folder(output_results_path, folder)
            titles = {}
            titles = amr.read_ascii_art(os.path.join(os.getcwd(), 'titles.txt'))
            print(titles['boosttimes'])
            for model in models:
                for k in range(8, 9):
                        print(f"Kmer Analysis: {k}")
                        print('\n\n'+type(model).__name__)
                        amr.process_model_results(output_results_path=output_results_path,
                                                folder=folder,
                                                model=model,
                                                kmer=k,
                                                antibiotic_dfs=antibiotic_dfs,
                                                input_kmer_path=input_kmer_path,
                                                entry=entry)
    elif entry == 'bond':
        for folder in tqdm(folders):
            amr.create_folder(output_results_path, folder)
            titles = {}
            titles = amr.read_ascii_art(os.path.join(os.getcwd(), 'titles.txt'))
            print(titles['bond'])
            for model in models:
                for k in range(8, 9):
                        print(f"Kmer Analysis: {k}")
                        print('\n\n'+type(model).__name__)
                        amr.process_model_results(output_results_path=output_results_path,
                                                folder=folder,
                                                model=model,
                                                kmer=k,
                                                antibiotic_dfs=antibiotic_dfs,
                                                input_kmer_path=input_kmer_path,
                                                entry=entry)
    elif entry == 'gexp':
        titles = {}
        titles = amr.read_ascii_art(os.path.join(os.getcwd(), 'titles.txt'))
        print(titles['models_exec_gexp'])
        amr.create_folder(output_results_path, entry)
        for model in models:
            print('\n\n'+type(model).__name__)
            amr.process_model_results(output_results_path=output_results_path,
                                      folder=entry,
                                      model=model,
                                      kmer=args.k_mer,
                                      antibiotic_dfs=antibiotic_dfs,
                                      input_kmer_path=input_kmer_path,
                                      entry=entry)

if __name__ == '__main__':
    main()
