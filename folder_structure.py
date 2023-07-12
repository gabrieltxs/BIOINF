import os
import amr_functions as amr
# script to generate default folder structure

# Set the root directory
root_dir = os.getcwd()
paths = ['lib\\kmer', 'lib\\gexp', 'lib\\bond', 'lib\\boost', 'lib\\target', 'lib\\genes_eda_data',
         'raw_data\\amr-phenotype', 'raw_data\\gexp-specialtygenes', 'raw_data\\ngs-specialtygenes','raw_data\\ngs-wgs',
         'scaffold',
         'results']
for path in paths:
    amr.create_folder(root_dir,path)
