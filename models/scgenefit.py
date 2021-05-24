import time

import numpy as np
from scGeneFit.functions import __sample, __select_constraints_centers, __select_constraints_pairwise, \
    __select_constraints_summarized, __lp_markers


def get_importance(data, labels, num_markers, method='centers', epsilon=1, sampling_rate=1, n_neighbors=3,
                   max_constraints=1000, redundancy=0.01, verbose=True):
    """marker selection algorithm
    data: Nxd numpy array with point coordinates, N: number of points, d: dimension
    labels: list with labels (N labels, one per point)
    num_markers: target number of markers to select. num_markers<d
    method: 'centers', 'pairwise', or 'pairwise_centers'
    epsilon: constraints will be of the form expr>Delta, where Delta is chosen to be epsilon times the norm of the smallest constraint (default 1)
    (This is the most important parameter in this problem, it determines the scale of the constraints,
    the rest the rest of the parameters only determine the size of the LP)
    sampling_rate: (if method=='pairwise' or 'pairwise_centers') selects constraints from a random sample of proportion sampling_rate (default 1)
    n_neighbors: (if method=='pairwise') chooses the constraints from n_neighbors nearest neighbors (default 3)
    max_constraints: maximum number of constraints to consider (default 1000)
    redundancy: (if method=='centers') in this case not all pairwise constraints are considered
    but just between centers of consecutive labels plus a random fraction of constraints given by redundancy
    if redundancy==1 all constraints between pairs of centers are considered """
    d = data.shape[1]
    t = time.time()
    samples, samples_labels, idx = __sample(data, labels, sampling_rate)

    if method == 'pairwise_centers':
        constraints, smallest_norm = __select_constraints_centers(
            data, labels, samples, samples_labels)
    elif method == 'pairwise':
        constraints, smallest_norm = __select_constraints_pairwise(
            data, labels, samples, samples_labels, n_neighbors)
    else:
        constraints, smallest_norm = __select_constraints_summarized(data, labels, redundancy)

    num_cons = constraints.shape[0]
    if num_cons > max_constraints:
        p = np.random.permutation(num_cons)[0:max_constraints]
        constraints = constraints[p, :]
    if verbose:
        print('Solving a linear program with {} variables and {} constraints'.format(
            constraints.shape[1], constraints.shape[0]))
    sol = __lp_markers(constraints, num_markers, smallest_norm * epsilon)
    if verbose:
        print('Time elapsed: {} seconds'.format(time.time() - t))
    return sol['x'][0:d]  # ndarray
