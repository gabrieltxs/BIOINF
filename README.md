# AMR prediction with LGBMClassifier models
This repository contains a Python script for predicting antimicrobial resistance (AMR) using the LGBMClassifier model. The script reads input datasets from a directory, applies feature extraction techniques to obtain k-mer features, trains and tests the models using cross-validation, and outputs the results in text files.


![PatricDB (5) vpd (1)-1](https://github.com/gabrieltxs/BIOINF/assets/43249674/3adce050-5492-4349-8504-7e7744466bfc)


## Getting Started
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites
This script requires the following Python libraries:

    pandas
    scikit-learn
    numpy
    tqdm
    lightgbm
    hyperopt
    joblib
    bayesian-optimization
    skopt

### Installing
Clone the repository to your local machine and install the required libraries:


```bash
  $ git clone https://github.com/username/repo.git
  $ cd repo
  $ pip install -r requirements.txt
```


### Usage
To use the script, execute the following command:

css
Copy code

```bash
  $ python main.py
```

## Code Structure
The main script consists of several sections:

    1 Import necessary libraries
    2 Set seed for reproducibility
    3 Define function to get list of models to evaluate
    4 Load list of selected samples
    5 Call function to get list of models
    6 Initialize KFold cross-validation
    7 Iterate over values of k to read the corresponding k-mer feature dataset
    8 Iterate over the models list
    9 Write results to text file

## Data Description
The input datasets are CSV files containing bacterial genomic sequences and their corresponding resistance profiles for selected antibiotics. The script reads these files from a directory and applies k-mer feature extraction techniques to obtain numerical feature vectors.

## Models
The script uses two models for AMR prediction: Random Forest and LGBMClassifier.

## Output
The script outputs the results of each model to a text file in the specified output directory. The results include accuracy, precision, recall, F1 score, and area under the ROC curve.

## Authors
Gabriel Sousa - gabrieltxs

## License
This project is licensed under the MIT License - see the LICENSE.md file for details.
[![MIT License](https://img.shields.io/badge/License-MIT-green.svg)](https://choosealicense.com/licenses/mit/)
