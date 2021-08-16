"""Detect communities using spatially-corrected modularity measure
then extract the backbone and do DCP on each community

All these should be run from Diss folder as python3.9 port_scripts/first_hierarchy.py
scale: γ
spatial null: gravity
null for DCP detection: ER (as Elliot et al.)
port network: weighted
"""

import os
from os.path import join
from tqdm import tqdm
import numpy as np
import scipy.sparse as sp
from scipy.io import savemat
import matlab.engine
from spatial_graphs import SpatialDiGraph
from spatial_graphs import util_funcs as f

# settings
weighting = "weighted"
year = "2019"  # ["2019", "2020", "2019/quarters/Q1"]
res: float = 0.85
spatial_null = "radiation"
constraint = "doubly"
# radiation 0.85 doesn't work

# load data
datadir = join("data", "matrices", year)
comdir = join("data", "port_results", year, weighting, "communities")
fmat = sp.load_npz(join(datadir, weighting, "cargo_ports.npz")).toarray()
dmat = np.load(join(datadir, "port_dmat.npy"))
g = SpatialDiGraph.from_numpy_array(fmat, dists=dmat)

B = getattr(f, f"modularity_{spatial_null}")(g, constraint=constraint, res=res)
B = (B + B.T) / 2.
eng = matlab.engine.start_matlab()
savemat("tempA.mat", {"array": g.fmat}, do_compression=False)
savemat("tempB.mat", {"array": B}, do_compression=False)
A = eng.load('tempA.mat')['array']
B = eng.load('tempB.mat')['array']
y, Q = eng.spectral23(A, B, nargout=2)
best = np.argmax(Q)
y, Q = [int(x) for x in y[best]], Q[0][best] / g.fmat.sum()
ncoms = len(set(y))
os.remove("tempA.mat")
os.remove("tempB.mat")

if spatial_null != "ng":
    np.savez(join(comdir, f"{spatial_null}_{constraint}_spectral_{res}"),
             partition=y, mod=Q, res=res, ncoms=ncoms)
else:
    np.savez(join(comdir, f"{spatial_null}_spectral_{res}"),
             partition=y, mod=Q, res=res, ncoms=ncoms)

print(f"found {ncoms} communities for resolution {res}.")
