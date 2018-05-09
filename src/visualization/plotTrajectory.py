import numpy as np
import random as random
import matplotlib.pylab as pylab
import matplotlib.pyplot as plt

# - Dimensions of interest.
dim1 = 1
dim2 = 0

params = {'legend.fontsize': 'x-large',
          'figure.figsize': (5, 4),
          'axes.labelsize': 10,
          'axes.titlesize': 'x-large',
          'xtick.labelsize': 10,
          'ytick.labelsize': 10}
pylab.rcParams.update(params)

fid = open('OUTPUT/inversion_1_trajecory.txt')
dummy = fid.read().strip().split()
fid.close()
dimensions = int(dummy[0])
iterations = int((dummy.__len__() - 2) / (dimensions + 1))

parameters = [[]]
for i in range(1, dimensions): parameters.append([])
dim1_model = []
dim2_model = []
misfit_model = []
for i in range(1, iterations + 1):
    dim1_model.append(float(dummy[2 + dim1 + (i - 1) * (dimensions + 1)]))
    dim2_model.append(float(dummy[2 + dim2 + (i - 1) * (dimensions + 1)]))
    misfit_model.append(dummy[2 + dimensions + (i - 1) * (dimensions + 1)])

plt.plot(dim1_model, dim2_model, 'k', linewidth=1, zorder=1)
plt.scatter(dim1_model, dim2_model, c=misfit_model, edgecolors='none', zorder=2, s=20)

plt.xlabel('q' + str(dim1 + 1))
plt.ylabel('q' + str(dim2 + 1))
plt.tight_layout()
plt.show()
# plt.savefig('OUTPUT/trajectory.pdf', format='pdf')
# plt.savefig('OUTPUT/trajectory.png')
