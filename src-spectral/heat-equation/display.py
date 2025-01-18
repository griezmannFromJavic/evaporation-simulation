import csv
import matplotlib.pyplot as plt
import numpy as np

# Function to read CSV and convert to 2D array
def read_csv(filename):
    with open(filename, 'r') as file:
        reader = csv.reader(file)
        data = [list(map(float, row)) for row in reader]
    return np.array(data)

# Read the CSV data
filename = 'solution.csv'
data = read_csv(filename)

# Plot the 2D result
plt.imshow(data, cmap='viridis', interpolation='nearest')
plt.colorbar(label='Value')
plt.title('2D Result from CSV')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.show()
