"""MDCP Benchmark development

N must be divisible by 8."""
import numpy as np
import matplotlib.pyplot as plt

p1 = 0.5
p2 = 0.5
rho = 1
N = 16

L = int(rho * N * (N - 1))
nodes = [x for x in range(N)]

# assign to binary communities
n = N // 2
part = np.zeros(N)
part[0:n] = 0
part[n:N] = 1

# idealised community block matrix
# matrix of {0,1} for communities
comm_mat = np.zeros([N, N])  # noqa
comm_mat[0:n, 0:n] = 1.
comm_mat[n:N, n:N] = 1.

# assign to cores and peripheries
n = N // 8
cp = np.zeros(N)
for i in range(8):
    cp[i * n: i * n + n] = i % 4

# idealised CP block matrix
# matrix of {0,1} for 2 CP structs
cp_mat = np.zeros([N, N])  # noqa
# put 1 inside block structures
M[:, n: 2 * n] = 1.
M[2 * n: 3 * n, :] = 1.
M[:, :n] = 0.
M[3 * n:, :] = 0.
