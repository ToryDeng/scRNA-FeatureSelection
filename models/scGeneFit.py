import numpy as np
from sklearn.neighbors import NearestNeighbors
from scipy.optimize import linprog


def select_constraints(data, labels, k):
    """
    Build constraints for linear programing.

    :param data: a N × d array (d dimension, N number of points)
    :param labels: a N × 1 array with labels from 1 to L
    :param k: the number of neighbors to select
    :return:
        Delta: a d × n array, where each column is of the form (v-w).^2 where v and w are nearest neighbors with different labels
        smallest: the smallest norm of vectors in Delta
    """
    n_total, d = data.shape
    smallest = np.inf
    classes = np.unique(labels)
    seperated = [data[np.where(labels == cell_class)] for cell_class in classes]  # 0, 1,..., L
    idx = 0

    Delta = np.zeros(shape=(k * n_total, d))
    for c in range(len(seperated)):
        seperated_copy = seperated.copy()
        del seperated_copy[c]
        current_class = seperated[c]
        other_classes = np.concatenate(seperated_copy)
        distances, indices = NearestNeighbors(n_neighbors=k, metric='euclidean').fit(other_classes).kneighbors(
            current_class)
        if smallest > distances.min():
            smallest = distances.min()
        for i in range(indices.shape[0]):
            for j in range(indices.shape[1]):
                nn_idx = indices[i, j]
                delta = current_class[i, :] - other_classes[nn_idx, :]
                Delta[idx, :] = delta
                idx += 1
    Delta, smallest = Delta * Delta, smallest * smallest
    return Delta.T, smallest


def sqz_hinge(Delta, epsilon, dim_target, lam=None):
    """
    Execute the linear programing and calculate importance for each gene.

    :param Delta: a d x N array of constraints, d is the dimension of the space and N is the number of constraints
    :param epsilon: a hinge parameter
    :param dim_target: number of markers to select
    :param lam: a N x d array stating a penalization for each constraint
        (typically all ones or inversely proportional to the number of on the number of points in the cluster)
    :return: M: a dim x 1 binary array. M(t)=1 implies that t is a selected marker
    """
    d, N = Delta.shape
    f = np.concatenate([np.zeros(shape=(d,)), np.ones(shape=(N,))])

    if lam is None:
        b = -epsilon * np.ones(shape=(N,))
    else:
        b = -epsilon * np.ones(shape=(N,)) / lam
    b = np.concatenate([b, np.ones((1,)) * dim_target], axis=0)

    A = np.concatenate([-Delta.T, -np.identity(N)], axis=1)
    bottom = np.concatenate([np.ones(shape=(1, d)), np.zeros(shape=(1, N))], axis=1)
    A = np.concatenate([A, bottom], axis=0)

    lb = np.zeros(shape=(d + N,))
    ub = np.concatenate([np.ones((d,)), np.inf * np.ones((N,))])
    bounds = [(left, right) for left, right in zip(lb, ub)]

    res = linprog(c=f, A_ub=A, b_ub=b, bounds=bounds, options={'tol': 1e-4})
    M = res['x'][:d]  # numpy.ndarray, just select all genes
    return M


def scGeneFit(data, labels, num_markers, opt=None):
    """
    The main function of scGeneFit.

    :param data: count matrix
    :param labels: cell types
    :param num_markers: number of markers to select
    :param opt: None
    :return: gene importance
    """
    constraints_neighbors = 4
    hinge_scale = 0.1

    D, s = select_constraints(data=data, labels=labels, k=constraints_neighbors)
    M = sqz_hinge(Delta=D, epsilon=hinge_scale * s, dim_target=num_markers)
    return M  # feature importance? hierarchy?
