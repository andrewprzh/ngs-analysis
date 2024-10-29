import sys

import matplotlib.pyplot as plt
import csv

import numpy

# Replace 'your_file.tsv' with the path to your file
file_path = sys.argv[1]

x_values = []
y_values = []

# Reading the TSV file
with open(file_path, 'r') as file:
    reader = csv.reader(file, delimiter='\t')
    next(reader)  # Skip the first line (header)

    for row in reader:
        try:
            # Read integers from the 2nd and 3rd columns
            x_values.append(int(row[1]))
            y_values.append(int(row[2]))
        except ValueError:
            print("Non-integer value encountered, skipping row:", row)

print(numpy.corrcoef(x_values, y_values)[0][1])

# Plotting the scatter plot
plt.figure(figsize=(8, 8))
plt.scatter(x_values, y_values, color='black', marker='o')
ax = plt.gca()
ax.set_xlim([0, 2000])
ax.set_ylim([0, 2000])
#ax.set_yscale('log')
#ax.set_xscale('log')
#plt.title('Scatter Plot of Values from TSV')
plt.xlabel('ONT')
plt.ylabel('Illumina')
plt.savefig(sys.argv[1] + '.png')