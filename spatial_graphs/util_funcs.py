"""Functions I use a lot in different scripts...needs tidying up etc."""

from spatial_graphs import SpatialGraph, SpatialDiGraph
from spatial_nets import GravityModel, RadiationModel, LocationsDataClass
import spatial_nets
import numpy as np


# utils
def symmetrise_digraph(g):
    fmat = (g.fmat + g.fmat.T) / 2
    dmat = g.dists
    g = SpatialGraph.from_numpy_array(fmat, dmat)
    return g

# modularities


def modularity_gravity(g, res=1, constraint=None,  **kwargs):
    locs = LocationsDataClass(g)
    if constraint == "unconstrained":
        grav = GravityModel(constraint=None)
        prediction = grav.fit_transform(locs)
        const = getattr(
            spatial_nets, f"{constraint.capitalize()}Model")()
    elif constraint in ["production", "attraction", "doubly"]:
        grav = GravityModel(constraint=constraint)
        prediction = grav.fit_transform(locs)
        const = getattr(
            spatial_nets, f"{constraint.capitalize()}Constrained")()
    else:
        raise ValueError("invalid constraint supplied")
    const_grav = const.fit_transform(locs, prediction)
    B = g.fmat - res * const_grav
    return B


def modularity_radiation(g, res=1, constraint=None, **kwargs):
    locs = LocationsDataClass(g)
    if constraint == "unconstrained":
        grav = RadiationModel(constraint=None)
        prediction = grav.fit_transform(locs)
        const = getattr(
            spatial_nets, f"{constraint.capitalize()}Model")()
    elif constraint in ["production", "attraction", "doubly"]:
        grav = RadiationModel(constraint=constraint)
        prediction = grav.fit_transform(locs)
        const = getattr(
            spatial_nets, f"{constraint.capitalize()}Constrained")()
    else:
        raise ValueError("invalid constraint supplied")
    const_grav = const.fit_transform(locs, prediction)
    B = g.fmat - res * const_grav
    return B


def modularity_ng(g, res=1, **kwargs):
    """"wrapper function for modularity_spa_GN"""
    O_vec = g.fmat.sum(axis=1)
    _, B, _ = modularity_spa_GN(g.fmat, g.dists, O_vec, 2)
    return B


def modularity_spa_GN(flow_mat, dist_mat, N, binsize, res=1):
    """
    Originally by Paul Expert.

    Arguments
    ---------
    flow_mat : adjacency matrix
    dist_mat : distance matrix between the nodes
    N : a measure of te importance of a node (by defaults its strength: Dist=sum(flow_mat,1) for example)
    binsize : size of the bins in the estimation of the deterrence function (has to be tuned)
    res : (my addition)
        allows resolution to be tuned

    Returns
    -------
    modularity_Spa
    modularity_GN
    deterrence_fct
    """

    nb_nodes, _ = flow_mat.shape
    nbins = int(np.ceil(dist_mat.max() / binsize)) + 1  #  I changed this

    deterrence_fct = np.zeros(nbins)
    norma_deterrence = np.zeros(nbins)

    matrix_distance = np.zeros((nb_nodes, nb_nodes), dtype=int)
    # null_model_GN = np.zeros_like(matrix_distance)
    null_model_Spa = np.zeros((nb_nodes, nb_nodes))

    # flow_mat = flow_mat + flow_mat.T  # symmetrised matrix (doesn't change the outcome of comm. detect: arXiv:0812.1770)
    degree = flow_mat.sum(axis=0)
    null_N = N[:, np.newaxis] * N
    matrix = flow_mat / null_N  # normalised adjacency matrix

    for i in range(nb_nodes):
        for j in range(nb_nodes):

            # convert distances in binsize's units
            dist = int(np.ceil(dist_mat[i, j] / binsize))
            matrix_distance[i, j] = dist

            # weighted average for the deterrence function
            num = matrix[i, j]
            deterrence_fct[dist] += num * N[i] * N[j]
            norma_deterrence[dist] += N[i] * N[j]

    # normalisation of the deterrence function
    for i in range(nbins):
        if norma_deterrence[i] > 0:
            deterrence_fct[i] /= norma_deterrence[i]

    # computation of the randomised correlations (preserving space), spatial null-model
    for i in range(nb_nodes):
        for j in range(nb_nodes):
            null_model_Spa[i, j] = deterrence_fct[matrix_distance[i, j]]

    # the modularity matrix for the spatial null-model
    prod = null_N * null_model_Spa
    modularity_Spa = flow_mat - res * prod * flow_mat.sum() / prod.sum()

    # the modularity matrix for the GN null-model
    null_model_GN = degree[:, np.newaxis] * degree / degree.sum()
    modularity_GN = flow_mat - res * null_model_GN

    return modularity_Spa, modularity_GN, deterrence_fct


# backbones
def get_gravity_backbone(g, constraint="unconstrained"):
    """Calculate positive and negative backbones of graph.

    TODO: comment each line explaining it
    Parameters:
    -----------
    g : spatial_graphs.SpatialGraph or spatial_graphs.SpatialDiGraph
    constraint : str
        one of ["unconstrained", "production", "attraction", "doubly"]

    References
    ----------
    ..[1] Rodrigo Leal Cervantes ... 
    """
    locs = LocationsDataClass(g)
    if constraint == "unconstrained":
        grav = GravityModel(constraint=None)
        prediction = grav.fit_transform(locs)
        constr = getattr(
            spatial_nets, f"{constraint.capitalize()}Model")()
    elif constraint in ["production", "attraction", "doubly"]:
        grav = GravityModel(constraint=constraint)
        prediction = grav.fit_transform(locs)
        constr = getattr(
            spatial_nets, f"{constraint.capitalize()}Constrained")()
    else:
        raise ValueError("invalid constraint supplied")
    _ = constr.fit_transform(locs, prediction)  # make it constrained
    pvals = constr.pvalues()
    pvals.set_significance()
    pos, neg = pvals.compute_backbone()  # returns matrices
    pos = pos.toarray() * 1
    neg = neg.toarray() * 1
    dists = g.dists
    pos_bb = SpatialDiGraph.from_numpy_array(fmat=pos, dists=dists)
    neg_bb = SpatialDiGraph.from_numpy_array(fmat=neg, dists=dists)
    #pos_bb = type(g).from_numpy_array(fmat=pos, dists=dists)
    #neg_bb = type(g).from_numpy_array(fmat=neg, dists=dists)
    return pos_bb, neg_bb


def get_radiation_backbone(g, constraint="unconstrained"):
    """Calculate positive and negative backbones of graph.

    TODO: comment each line explaining it
    Parameters:
    -----------
    g : spatial_graphs.SpatialGraph or spatial_graphs.SpatialDiGraph
    constraint : str
        one of ["unconstrained", "production", "attraction", "doubly"]

    References
    ----------
    ..[1] Rodrigo Leal Cervantes ... 
    """
    locs = LocationsDataClass(g)
    if constraint == "unconstrained":
        rad = GravityModel(constraint=None)
        prediction = rad.fit_transform(locs)
        constr = getattr(
            spatial_nets, f"{constraint.capitalize()}Model")()
    elif constraint in ["production", "attraction", "doubly"]:
        rad = RadiationModel(constraint=constraint)
        prediction = rad.fit_transform(locs)
        constr = getattr(
            spatial_nets, f"{constraint.capitalize()}Constrained")()
    else:
        raise ValueError("invalid constraint supplied")
    _ = constr.fit_transform(locs, prediction)  # TODO: explain this line
    pvals = constr.pvalues()
    pvals.set_significance()
    pos, neg = pvals.compute_backbone()  # returns matrices
    pos = pos.toarray() * 1
    neg = neg.toarray() * 1
    dists = g.dists
    pos_bb = type(g).from_numpy_array(fmat=pos, dists=dists)
    neg_bb = type(g).from_numpy_array(fmat=neg, dists=dists)
    return pos_bb, neg_bb
