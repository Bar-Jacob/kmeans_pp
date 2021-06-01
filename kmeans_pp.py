import pandas as pd
import numpy as np
import sys
import random
#import mykmeanssp


k = sys.argv[1]


def argu_check(k, max_iter=200):
    if(('.' in k) or int(k) < 0):
        print("Invalid Input")
        assert()
    if(max_iter != 200):
        if(('.' in max_iter) or int(max_iter) < 0):
            print("Invalid Input")
            assert()
    return int(max_iter)


try:
    max_iter = sys.argv[2]
    max_iter = argu_check(k, max_iter)
    input1 = sys.argv[3]
    input2 = sys.argv[4]
except IndexError:
    max_iter = argu_check(k)
    input1 = sys.argv[2]
    input2 = sys.argv[3]

k = int(sys.argv[1])


points1 = pd.read_csv(input1)
points2 = pd.read_csv(input2)
points1.columns = ['c1', 'c2', 'c3']
points2.columns = ['c1', 'c2', 'c3']
merged_points = pd.merge(points1, points2, on='c1')
merged_points = merged_points.set_index("c1")
merged_points = merged_points.sort_index(ascending=True)
print(merged_points)

centroids = [0.0 for i in range(k)]
distances = [0.0 for i in range(len(merged_points))]
probabilities = [0.0 for i in range(len(merged_points))]

merged_points = merged_points.values.tolist()
random_index = random.randint(0, 97)


def probabilities_calc():
    sum = 0.0
    for distance in distances:
        sum = sum + distance
    for i in range(len(distances)):
        probabilities[i] = distances[i]/sum


def min_distances(Z):
    cnt = 0
    for vector in merged_points:
        min = -1
        for i in range(Z):
            sum = 0.0
            for xi in range(len(vector)):
                sum += pow((vector[xi]-centroids[i][xi]), 2)
            if(sum < min or min == -1):
                min = sum
        distances[cnt] = min
        cnt += 1


def KMeansPP():
    Z = 1
    centroids[0] = merged_points[random_index]
    while(Z != k):
        min_distances(Z)
        probabilities_calc()
        random_vector = random.choices(merged_points, probabilities)
        centroids[Z] = random_vector[0]
        Z = Z + 1


KMeansPP()
print(centroids)
print(merged_points.index(centroids[0]))

#print(mykmeanssp.hello(1, 3.2))
