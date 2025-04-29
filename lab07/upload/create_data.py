from sklearn.datasets import make_blobs
from sklearn.preprocessing import StandardScaler
import numpy as np

centers = [[1, 1], [-1, -1], [1, -1], [-1.5, -1.5], [-2, 2], [1, 3]]
X, labels_true = make_blobs(
    n_samples=27*1024, centers=centers, cluster_std=0.25, random_state=0
)

X = StandardScaler().fit_transform(X)

np.savetxt("data", X)
