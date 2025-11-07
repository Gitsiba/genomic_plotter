import matplotlib.pyplot as plt
import numpy as np

# Simulate mutation impact scores across 5 domains
domains = ["N-term", "DBD", "Linker", "CTD", "Tail"]
impact_scores = np.random.rand(len(domains)) * 100

# Plot
plt.figure(figsize=(8, 5))
bars = plt.bar(domains, impact_scores, color="salmon")
plt.title("Mutation Impact Across Protein Domains", fontsize=14)
plt.ylabel("Impact Score", fontsize=12)
plt.tight_layout()
plt.savefig("../plots/mutation_impact_plot.png", dpi=300)
plt.show()
