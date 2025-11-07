import matplotlib.pyplot as plt
import numpy as np

# Simulated chromosome names and mutation counts
chromosomes = [f"chr{i}" for i in range(1, 13)]
mutation_counts = np.random.randint(50, 300, size=len(chromosomes))

# Normalize mutation counts for radial plotting
max_count = max(mutation_counts)
normalized = mutation_counts / max_count

# Plot
fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))
theta = np.linspace(0.0, 2 * np.pi, len(chromosomes), endpoint=False)
bars = ax.bar(theta, normalized, width=0.4, color='mediumorchid', alpha=0.8)

# Add labels
for i, bar in enumerate(bars):
    rotation = np.degrees(theta[i])
    ax.text(theta[i], normalized[i] + 0.05, chromosomes[i],
            ha='center', va='bottom', rotation=rotation,
            rotation_mode='anchor', fontsize=10)

ax.set_title("Circos-style Mutation Density Plot", fontsize=14)
ax.set_yticklabels([])
ax.set_xticklabels([])
plt.tight_layout()
plt.savefig("../plots/circos_mutation_plot.png", dpi=300)
plt.show()
