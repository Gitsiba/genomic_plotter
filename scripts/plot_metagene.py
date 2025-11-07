import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.ndimage import gaussian_filter1d

# Simulate synthetic expression data: 50 genes × 100 bins
np.random.seed(42)
num_genes = 50
num_bins = 100
expression_data = np.random.rand(num_genes, num_bins) * np.linspace(0.5, 1.5, num_bins)

# Compute average signal across genes
mean_signal = expression_data.mean(axis=0)

# Apply Gaussian smoothing
smoothed_signal = gaussian_filter1d(mean_signal, sigma=2)

# Plotting
sns.set(style="whitegrid")
plt.figure(figsize=(10, 6))
plt.plot(range(1, num_bins + 1), smoothed_signal, color="darkblue", linewidth=2)
plt.title("Metagene Plot: Average Signal Across Gene Bodies", fontsize=14)
plt.xlabel("Normalized Gene Body (TSS → TES)", fontsize=12)
plt.ylabel("Average Signal Intensity", fontsize=12)
plt.tight_layout()
plt.savefig("../plots/metagene_plot.png", dpi=300)
plt.show()
