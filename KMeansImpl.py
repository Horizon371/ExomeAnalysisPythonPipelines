import numpy as np
from numpy.linalg import norm


def find_closest_cluster(distances):
    return np.argmin(distances, axis=1)


class KMeansImpl:

    def __init__(self, k, max_iterations=100):
        self.k = k
        self.max_iterations = max_iterations
        self.labels = None
        self.centroids = None

    def initialize_centroids(self, X):
        random_idx = np.random.permutation(X.shape[0])
        centroids = X[random_idx[:self.k]]
        return centroids

    def compute_centroids(self, X, labels):
        centroids = np.zeros((self.k, X.shape[1]))
        for k in range(self.k):
            centroids[k, :] = np.mean(X[labels == k, :], axis=0)
        return centroids

    def compute_distances(self, X, centroids):
        distances = np.zeros((X.shape[0], self.k))
        for k in range(self.k):
            row_norm = norm(X - centroids[k, :], axis=1)
            distances[:, k] = row_norm
        return distances

    def fit(self, X):
        self.centroids = self.initialize_centroids(X)
        for i in range(self.max_iterations):
            old_centroids = self.centroids
            distances = self.compute_distances(X, old_centroids)
            self.labels = find_closest_cluster(distances)
            self.centroids = self.compute_centroids(X, self.labels)
            if np.all(old_centroids == self.centroids):
                break

    def predict(self, X):
        distances = self.compute_distances(X, self.centroids)
        return find_closest_cluster(distances)
