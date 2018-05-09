import matplotlib.pylab as pylab
import numpy as np
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

params = {'legend.fontsize': 'x-large',
          'figure.figsize': (12, 4),
          'axes.labelsize': 16,
          'axes.titlesize': 'x-large',
          'xtick.labelsize': 16,
          'ytick.labelsize': 16}
pylab.rcParams.update(params)

# If there's other stuff on the line, like misfit values
trailingElements = 1
nbi = 50
fid = open(sys.argv[1])
dummy = fid.read().strip().split()
fid.close()
numParameters = int(dummy[0])
numSamples = int(dummy.__len__() - 2) / (numParameters + 1)
parameters = []
for i in range(0, numParameters):
    parameters.append([])

for i in range(nbi, numSamples):
    for index, parameter in enumerate(parameters):
        parameter.append(float(dummy[2 + index + (i * (numParameters + trailingElements))]))

means = []
for index, parameter in enumerate(parameters):
    means.append(np.mean(parameter))

variances = []
for index, parameter in enumerate(parameters):
    variances.append(np.var(parameter))

maxVar = np.max(np.abs(variances))
f, axarr = plt.subplots(1, 2)
axarr[0].scatter(np.linspace(1,np.array(means).size,np.array(means).size),means)
axarr[0].set(xlabel='parameter #', ylabel='mean')
axarr[1].scatter(np.linspace(1,np.array(variances).size,np.array(variances).size),variances)
axarr[1].set(xlabel='parameter #', ylabel='standard deviation')

plt.show()