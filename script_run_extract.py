import subprocess

def run_kmer_extraction(file_path, output_path, k_value, gene_type, core_num):
    command = f"python A_kmerExtract.py -f main -k {k_value} -p {file_path} -o {output_path} -t {gene_type} -c {core_num}"
    print(f"Running command: {command}")
    subprocess.call(command, shell=True)
    print(f"Command completed: {command}")

# Modify the arguments as needed
#    "scaffold_genes/antibiotic_resistance",
#    "scaffold_genes/drug_target",
file_paths = [
    "scaffold_genes/antibiotic_resistance",
    "scaffold_genes/drug_target",
    "scaffold_genes/transporter",
    "scaffold_genes/virulence_factor"
]
output_path = "lib/kmer"
k_value = 8
core_num = 16
#"amr", "dt", 
gene_types = ["amr", "dt", "tpt", "vf"]

total_commands = len(file_paths)
completed_commands = 0

for i, file_path in enumerate(file_paths):
    gene_type = gene_types[i]
    run_kmer_extraction(file_path, output_path, k_value, gene_type, core_num)
    completed_commands += 1
    progress = (completed_commands / total_commands) * 100
    print(f"Progress: {progress:.2f}%")

command2 = f"python E_MLmodelsExec.py -k 5 -pk lib\kmer -pf lib\\processed -o results"
subprocess.call(command2, shell=True)

print("All commands executed.")