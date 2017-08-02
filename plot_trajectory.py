import numpy as np
import random as random
import matplotlib.pylab as pylab
import matplotlib.pyplot as plt

# ============================================================
# - Setup.
# ============================================================

# - Dimensions of interest.
dim1 = 0
dim2 = 4

params = {'legend.fontsize': 'x-large',
          'figure.figsize': (4, 4),
          'axes.labelsize': 10,
          'axes.titlesize': 'x-large',
          'xtick.labelsize': 10,
          'ytick.labelsize': 10}
pylab.rcParams.update(params)

# ============================================================
# - Read samples and plot trajectory.
# ============================================================

fid = open('OUTPUT/trajectory.txt')
dummy = fid.read().strip().split()
fid.close()
dimensions = int(dummy[0])
iterations = (dummy.__len__() - 2) / (dimensions + 1)


parameters = [[]]
for i in range(1, dimensions): parameters.append([])
dim1_model = []
dim2_model = []
misfit_model = []
for i in range(1, iterations + 1):
    dim1_model.append(float(dummy[2 + dim1 + (i - 1) * (dimensions + 1)]))
    dim2_model.append(float(dummy[2 + dim2 + (i - 1) * (dimensions + 1)]))
    # dim1_model[i-1] = 1/dim1_model[i-1]
    # dim2_model[i-1] = 1/dim2_model[i-1]
    misfit_model.append(dummy[2 + dimensions + (i - 1) * (dimensions + 1)])

plt.plot(dim1_model, dim2_model, 'k', linewidth=1, zorder=1)
plt.scatter(dim1_model, dim2_model, c=misfit_model, edgecolors='none', zorder=2, s=20)

# # Plotting gradient
# fid = open('OUTPUT/gradient.txt')
# dummy = fid.read().strip().split("\n")
# fid.close()
#
# # parse data
# q1 = []
# q2 = []
# dq1 = []
# dq2 = []
# for line in dummy:
#     splitty = line.split()
#     q1.append(float(splitty[0]))
#     q2.append(float(splitty[1]))
#     dq1.append(float(splitty[2]))
#     dq2.append(float(splitty[3]))
# Q = plt.quiver(q1, q2, dq1, dq2)

plt.xlabel('q' + str(dim1 + 1))
plt.ylabel('q' + str(dim2 + 1))
# plt.gca().set_aspect('equal', adjustable='box')
# plt.xlim([0.5, 2.0])
# plt.ylim([2.65, 3.4])
plt.savefig('OUTPUT/trajectory.pdf', format='pdf')
# plt.savefig('OUTPUT/trajectory.png')
