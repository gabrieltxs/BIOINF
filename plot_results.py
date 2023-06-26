import os
import re
import matplotlib.pyplot as plt
import numpy as np

base_path = os.path.join(os.getcwd(), 'results')  # base folder path

folders = ['amr', 'dt', 'tpt', 'vf', 'wgs', 'gpa']
pattern = r'(\d+\.\d+);(\d+\.\d+);'
k = 8
# Create a figure and subplots for each folder
fig, axs = plt.subplots(len(folders), 1, figsize=(8, 6 * len(folders)), sharex=True)

# Traverse the folders
for i, folder in enumerate(folders):
    folder_path = os.path.join(base_path, folder)

    validation_scores = []
    test_scores = []

    # Traverse the files in the folder
    for file_name in os.listdir(folder_path):
        file_path = os.path.join(folder_path, file_name)
        if re.search(str(k), file_name):

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

    # Plot the scores in the respective subplot
    ax = axs[i]
    x = np.arange(len(validation_scores))
    width = 0.4

    ax.bar(x - width/2, validation_scores, width, align='center', alpha=0.8, label='Validation F1-score')
    ax.bar(x + width/2, test_scores, width, align='center', alpha=0.8, label='Test F1-score')

    ax.set_ylim(0.4, 1)
    ax.set_ylabel('F1-score')
    ax.set_title(folder.capitalize())
    category = ['ceftazidime', 'ciprofloxacin', 'meropenem', 'tobramycin']
    ax.set_xticks(x)
    ax.set_xticklabels(category)
    ax.legend()

    # Add value labels to each bar
    for j, v in enumerate(validation_scores):
        ax.text(x[j] - width/2, v + 0.01, str(round(v, 3)), ha='center', va='bottom')
    for j, v in enumerate(test_scores):
        ax.text(x[j] + width/2, v + 0.01, str(round(v, 3)), ha='center', va='bottom')

# Adjust spacing between subplots
plt.tight_layout()

# Save the plot as an image file
output_file = os.path.join(os.getcwd(), f"scores_plot{str(k)}.png")
plt.savefig(output_file)
