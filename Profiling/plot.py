import matplotlib.pyplot as plt
import numpy as np

# Execution times for the two versions
a = "15.796 14.6274 15.789 16.5872 15.293 15.9239 15.67 16.0072 15.4953 14.4392"
a = [float(x) for x in a.split()]

b = "15.0694 14.3055 14.4349 14.1407 14.254 15.384 13.99 14.8507 15.0578 13.8498"
b = [float(x) for x in b.split()]

# Create a box plot
data = [a, b]
plt.boxplot(data, labels=['FWaveVecSolver', 'FWaveVecSolver with fma'], showfliers=None)

# Add title and labels
plt.ylabel('Execution Time (s)')

# Show the plot
plt.grid()
plt.show()