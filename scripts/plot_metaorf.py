import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Simulate translation signal across 3 ORFs
np.random.seed(1)
orf_lengths = [150, 300, 120]
orf_signals = [np.random.rand(l) * np.linspace(1, 2, l) for l in orf_lengths]

# Plot
sns.set(style="whitegrid")
plt.figure(figsize=(10, 6))
for i, signal in enumerate(orf_signals):
    plt.plot(range(1, len(signal) + 1), signal, label=f"ORF {i+1}")
plt.title("metaORF Plot: Translation Signal Across ORFs", fontsize=14)
plt.xlabel("Position in ORF", fontsize=12)
plt.ylabel("Signal Intensity", fontsize=12)
plt.legend()
plt.tight_layout()
plt.savefig("../plots/metaorf_plot.png", dpi=300)
plt.show()
