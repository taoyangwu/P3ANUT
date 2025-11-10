


import numpy as np
import matplotlib.pyplot as plt

data = np.load("testMerge.npy", allow_pickle=True)
data = np.unique(data, return_counts=True)
counts = data[1] / np.sum(data[1])  # Normalize to frequency
plt.bar(data[0], counts, width=1.0, alpha=0.7)
plt.xlabel('Number of Insert/Delete Operations')
plt.ylabel('Frequency (log scale)')
plt.title('Distribution of Insert/Delete Operations in Merged Reads')
plt.grid(True)
plt.show()
