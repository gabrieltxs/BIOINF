import os
import amr_functions as amr
# script to generate default folder structure

# Set the root directory
root_dir = os.getcwd()

# Create the directories if they don't already exist
if not os.path.exists(os.path.join(root_dir, "lib", "kmer")):
    os.makedirs(os.path.join(root_dir, "lib", "kmer"))
    print("Folder 'kmer' created successfully in 'lib'")
else:
    print("Folder 'kmer' already exists in 'lib'")

if not os.path.exists(os.path.join(root_dir, "lib", "target")):
    os.makedirs(os.path.join(root_dir, "lib", "target"))
    print("Folder 'target' created successfully in 'lib'")
else:
    print("Folder 'target' already exists in 'lib'")

if not os.path.exists(os.path.join(root_dir, "results")):
    os.makedirs(os.path.join(root_dir, "results"))
    print("Folder 'results' created successfully")
else:
    print("Folder 'results' already exists")

if not os.path.exists(os.path.join(root_dir, "scaffold")):
    os.makedirs(os.path.join(root_dir, "scaffold"))
    print("Folder 'scaffold' created successfully")
else:
    print("Folder 'scaffold' already exists")

amr.create_folder(root_dir, 'raw_data')