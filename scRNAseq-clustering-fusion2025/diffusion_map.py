import numpy as np
import pandas as pd
from pydiffmap import diffusion_map as dmap

# Step 1: Prepare or load your consensus matrix
# For example, create a dummy consensus matrix (replace this with your actual matrix)
# Ensure it is symmetric and represents similarity between clusters
cell_pair_counts = pd.read_csv("/mnt/pan/SOM_CCCC_JGS25/shultesp/scRNAseq-clustering-fusion2025/cell_pair_counts.csv", index_col=0)
# Convert to a distance matrix (1-normalized co-occurrence)
cooccurrence_matrix_normalized = cell_pair_counts.div(1000) #divide by number of iterations
distances = 1 - cooccurrence_matrix_normalized
distances
consensus_matrix=distances

# Step 2: Initialize and run Diffusion Map
dm = dmap.DiffusionMap(n_evecs=2)  # 'p' for the probabilistic diffusion map
evecs = dm.fit_transform(consensus_matrix)

# Step 3: Convert eigenvectors to a usable format
# You can convert the output to a DataFrame if needed
import pandas as pd

evecs_df = pd.DataFrame(evecs, columns=['DiffusionMap1', 'DiffusionMap2'])
print(evecs_df)

# Step 4: (Optional) Visualize the results
import matplotlib.pyplot as plt

plt.scatter(evecs_df['DiffusionMap1'], evecs_df['DiffusionMap2'])
plt.title('Diffusion Map of Consensus Matrix from Bootstrapping n=1000 times')
plt.xlabel('Diffusion Map 1')
plt.ylabel('Diffusion Map 2')
plt.savefig("consensus_matrix_diffusion.png")
