import pandas as pd
import matplotlib.pyplot as plt

# Path to your CSV file
csv_file_path = r'C:\BIOINF\BIOINF\test_scores.csv'

# Read the CSV file using pandas
df = pd.read_csv(csv_file_path, delimiter=';')

# Extract the column names excluding the first column
columns = df.columns[1:]

# Set the figure size to be four times wider
plt.figure(figsize=(20, 10))

# Iterate over each column and plot the values as scatter points
for i, column in enumerate(columns):
    color = plt.cm.tab10(i)  # Get a color from the tab10 colormap for each column
    plt.scatter(df['Folder'], df[column], label=column, color=color)

    # Connect the scatter points to the x-axis with dotted lines of the same color
    for x, y in zip(df['Folder'], df[column]):
        plt.plot([x, x], [0, y], linestyle='dotted', color=color)

    # Find the maximum value for the column
    max_value = df[column].max()
    
    # Draw a horizontal line at the maximum value with the same color
    plt.axhline(y=max_value, color=color, linestyle='--')

# Set plot labels and title
plt.xlabel('Boost Combination')
plt.ylabel('Values')
plt.title('F1-Score')

# Rotate x-axis labels by 45 degrees
plt.xticks(rotation=45)

# Set the y-axis limits
plt.ylim(0.5, 1)

# Add a legend
plt.legend()

# Adjust the subplot parameters to expand downwards
plt.subplots_adjust(bottom=0.2)

# Save the plot in the current folder
plt.savefig('lib\genes_eda_data\scatter_plot.png')


