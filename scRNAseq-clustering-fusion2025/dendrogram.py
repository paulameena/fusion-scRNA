import scanpy as sc
import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score
from scipy.cluster.hierarchy import linkage, dendrogram
import itertools
import matplotlib.pyplot as plt

# Hierarchical clustering using ward linkage (can experiment with other methods)
#linkage_matrix = linkage(distances, method='ward')
linkage_matrix = pd.read_csv("/mnt/pan/SOM_CCCC_JGS25/shultesp/scRNAseq-clustering-fusion2025/linkage.csv")
linkage_matrix
#linkage_matrix.index = distances.index


#The silhouette score is not suitable as a metric for hierarchical clustering, since it's a metric to asses cluster quality
#You could calculate the silhouette score on the final clusters after hierarchical clustering

#You can determine the number of clusters by looking at the dendrogram (see below)
#or by using some other metric (e.g. calculating the silhouette score for different numbers of clusters and selecting the one with the best score)
dendrogram(linkage_matrix)
plt.savefig("/mnt/pan/SOM_CCCC_JGS25/shultesp/scRNAseq-clustering-fusion2025/dendrogram.png")
