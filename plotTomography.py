import matplotlib.pylab as pylab
import numpy as np
import matplotlib.pyplot as plt
import math

params = {'legend.fontsize': 'x-large',
          'figure.figsize': (6, 6),
          'axes.labelsize': 16,
          'axes.titlesize': 'x-large',
          'xtick.labelsize': 16,
          'ytick.labelsize': 16}
pylab.rcParams.update(params)
#
# If there's other stuff on the line, like misfit values
trailingElements = 1
nbi = 30
fid = open('OUTPUT/samples1.txt')
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
plt.imshow(np.transpose(np.reshape(variances, (11,11))), cmap=plt.get_cmap('seismic'), interpolation='none',
           extent=[-0.5, 10 + 0.5, 10 + 0.5, -0.5], vmin=0, vmax=maxVar)
cbar = plt.colorbar(ticks=np.arange(0, maxVar, maxVar/10))
# plt.gca().invert_yaxis()
plt.savefig("OUTPUT/tomographyVariances.pdf", format='pdf')
plt.close()

plt.imshow(1 / np.transpose(np.reshape(means, (11,11))), cmap=plt.get_cmap('binary'), interpolation='none',
           extent=[-0.5, 10 + 0.5, -0.5, 10 + 0.5], vmin=0000, vmax=3001)
cbar = plt.colorbar(ticks=np.arange(0, 3001, 500))
plt.gca().invert_yaxis()
cbar.set_label('Speed of sound [m/s]')
plt.savefig("OUTPUT/tomographyMeans.pdf", format='pdf')
plt.close()
