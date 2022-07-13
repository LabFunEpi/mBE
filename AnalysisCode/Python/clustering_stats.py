from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score, calinski_harabasz_score
import numpy as np
import sys

if len(sys.argv) <= 1:
    sys.exit("Usage: python clustering_stats.py <infile>")

infile = sys.argv[1]

X = np.loadtxt(infile, delimiter='\t')

for i in range(2,12):
    n_clusters = i
    clusterer = KMeans(n_clusters=n_clusters, random_state=10)
    cluster_labels = clusterer.fit_predict(X)
    silhouette_avg = silhouette_score(X, cluster_labels)
    ch = calinski_harabasz_score(X, cluster_labels)

    print(str(i) + "\t" + str(silhouette_avg) + "\t" + str(ch))
