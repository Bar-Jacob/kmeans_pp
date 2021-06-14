import pandas as pd
import numpy as np
import sys
import random
import mykmeanssp

k = sys.argv[1]


def argu_check(k, max_iter=300):
    if(('.' in k) or int(k) < 0):
        print("Invalid Input")
        assert()
    if(max_iter != 300):
        if(('.' in max_iter) or int(max_iter) < 0):
            print("Invalid Input")
            assert()
    return int(max_iter)


try:
    check_num_of_args = sys.argv[4]
    max_iter = sys.argv[2]
    max_iter = argu_check(k, max_iter)
    input1 = sys.argv[3]
    input2 = sys.argv[4]
except IndexError:
    max_iter = argu_check(k)
    input1 = sys.argv[2]
    input2 = sys.argv[3]

k = int(k)

f1 = open(input1, 'r')
f2 = open(input2, 'r')
len_f1 = len(f1.readline().split(","))
len_f2 = len(f2.readline().split(","))
f1.close
f2.close

points1 = pd.read_csv(input1, names=["c" + str(i) for i in range(len_f1)])
points2 = pd.read_csv(input2, names=["c" + str(i) for i in range(len_f2)])

points1 = pd.DataFrame(points1)
points2 = pd.DataFrame(points2)

merged_points = points1.merge(points2, on='c0')
merged_points = merged_points.sort_values(by=["c0"])

distances = [-1.0 for i in range(len(merged_points))]
probabilities = [0.0 for i in range(len(merged_points))]

merged_points = merged_points.set_index('c0')
merged_points_numpy = merged_points.to_numpy()

centroids = np.array(
    [[0.0 for i in range(len(merged_points_numpy[0]))] for i in range(k)])
# creating a list of the locations of the centorinds in the merged_points_numpy list
centroids_locations = [0 for i in range(k)]

if(k >= len(merged_points_numpy)):
    print("There are more clusters than points")
    assert()


def probabilities_calc():
    sum = 0.0
    for distance in distances:
        sum += distance
    for i in range(len(distances)):
        probabilities[i] = distances[i]/sum


def min_distances(Z):
    for index in range(len(merged_points_numpy)):
        curr_distance = pow(np.linalg.norm(
            (merged_points_numpy[index]-centroids[Z-1])), 2)
        if(curr_distance < distances[index] or distances[index] == -1.0):
            distances[index] = curr_distance


def KMeansPP():
    cnt = 1
    np.random.seed(0)

    random_index = (int)(np.random.choice(merged_points.index))
    centroids_locations[0] = random_index
    centroids[0] = np.ndarray.copy(merged_points_numpy[random_index])
    Z = 1
    while(Z != k):
        min_distances(Z)
        probabilities_calc()
        random_index = (int)(np.random.choice(
            merged_points.index, p=probabilities))
        centroids_locations[cnt] = random_index
        centroids[Z] = np.ndarray.copy(merged_points_numpy[random_index])
        Z += 1
        cnt += 1


KMeansPP()
# setting an argument for C
dimension = len(centroids[0])
data_points_size = len(merged_points_numpy)

# creating a 1D list of the data points to give as an argument to C
data_points_p = []
for vector in merged_points_numpy:
    for i in range(dimension):
        data_points_p.append(vector[i])

final_centroids = np.array(mykmeanssp.fit(
    k, max_iter, dimension, data_points_size, centroids_locations, data_points_p))

print(*centroids_locations, sep=",")
for centroid in final_centroids:
    print(*[np.round(num, decimals=4) for num in centroid], sep=",")
