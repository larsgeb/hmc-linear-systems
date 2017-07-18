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
dimensions = int(dummy[0])
iterations = (dummy.__len__()-2)/(dimensions+1)

# Range from 1 to number of parameters
dim1 = 0
dim2 = 1

parameters = [[]]
for i in range(1, dimensions): parameters.append([])
dim1_model = []
dim2_model = []
misfit_model = []
for i in range(1, iterations + 1):
    dim1_model.append(dummy[2 + dim1 + (i - 1) * (dimensions+1)])
    dim2_model.append(dummy[2 + dim2 + (i - 1) * (dimensions+1)])
    misfit_model.append(dummy[2 + dimensions + (i - 1) * (dimensions+1)])


# plt.plot(dim1_model, dim2_model, 'k', linewidth=0.05)
plt.scatter(dim1_model, dim2_model,c=misfit_model)
plt.savefig('OUTPUT/trajectory.png')
