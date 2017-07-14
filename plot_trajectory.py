import numpy as np
import random as random
import matplotlib.pylab as pylab
import matplotlib.pyplot as plt

# ============================================================
# - Setup.
# ============================================================

# - Dimensions of interest.
dim_1 = 1
dim_2 = 2
# - Incremental displacement for duplicate points.
epsilon_1 = 0.0003
epsilon_2 = 0.0003

params = {'legend.fontsize': 'x-large',
          'figure.figsize': (8, 8),
          'axes.labelsize': 20,
          'axes.titlesize': 'x-large',
          'xtick.labelsize': 20,
          'ytick.labelsize': 20}
pylab.rcParams.update(params)

# ============================================================
# - Read samples and plot trajectory.
# ============================================================

fid = open('OUTPUT/trajectory.txt')
dummy = fid.read().strip().split()
fid.close()
print
dimensions = int(dummy[0])
iterations = int(dummy[1]) - 30

# Range from 1 to number of parameters
dim1 = 9
dim2 = 2

parameters = [[]]
for i in range(1, dimensions): parameters.append([])
dim1_model = []
dim2_model = []
dim1_momentum = []
dim2_momentum = []
for i in range(1, iterations + 1):
    dim1_model.append(dummy[1 + dim1 + 2 * (i - 1) * dimensions])
    dim1_momentum.append(dummy[1 + dim1 + 2 * (i) * dimensions])
    dim2_model.append(dummy[1 + dim2 + 2 * (i - 1) * dimensions])
    dim2_momentum.append(dummy[1 + dim2 + 2 * (i) * dimensions])

plt.plot(dim1_model, dim2_model, 'k', linewidth=0.05)
plt.plot(dim1_model, dim2_model, 'ro', linewidth=0.05, markersize=0.5)
plt.savefig('OUTPUT/trajectory.png')
