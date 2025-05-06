import pylab as plt
import numpy as np

plt.figure()
points = plt.loadtxt("clustered")
cluster_index_column = 2
clusters = np.unique(points[:, cluster_index_column])
print(clusters)
for c in clusters:
    points_in_cluster = points[np.where(
        points[:, cluster_index_column] == c)[0]]
    plt.scatter(points_in_cluster[:, 0], points_in_cluster[:, 1], label=c)

plt.show()
