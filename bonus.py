from sklearn import datasets
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import numpy as np
import matplotlib as mpl
mpl.use('agg')

# inspiration was taken from https://stackoverflow.com/questions/41540751/sklearn-kmeans-equivalent-of-elbow-method
data_set = datasets.load_iris()
data_frame = data_set.data

inertia_list = []
for k in range(1, 10):
    kmeans = KMeans(n_clusters=k, init='k-means++',
                    random_state=0)
    kmeans.fit(data_frame)
    inertia_list.append(kmeans.inertia_)

inertia_plot = plt.figure()
plt.plot(range(1, 10), inertia_list)
plt.title('Elbow method for selection of optimal "k" clusters')
plt.xlabel('K')
plt.ylabel('Average Dispersion')

plt.annotate("Elbow point", (2.3, inertia_list[1]), xytext=(
    3, 170), arrowprops=dict(facecolor='pink', shrink=0.05, linestyle="--"))
plt.scatter(2, inertia_list[1], s=500,
            facecolors='none', edgecolors='k', linestyle="--")
plt.savefig("elbow.png")
