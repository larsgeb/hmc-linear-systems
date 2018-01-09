import numpy as np
import sys

import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
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
nbi = 30
fid = open(sys.argv[1])
dummy = fid.read().strip().split()
fid.close()
numParameters = int(dummy[0])
numSamples = int((dummy.__len__() - 2) / (numParameters + 1))
parameters = []
for i in range(0, numParameters):
    parameters.append([])

for i in range(nbi, numSamples):
    for index, parameter in enumerate(parameters):
        parameter.append(float(dummy[2 + index + (i * (numParameters + trailingElements))]))

means = []
variances = []
for index, parameter in enumerate(parameters):
    means.append(np.mean(parameter))
    variances.append(np.var(parameter))

# for index, parameter in enumerate(parameters):

maxVar = np.max(np.abs(variances))
f, axarr = plt.subplots(1, 2)
image0 = axarr[0].imshow(
    np.flipud(np.transpose(np.reshape(means, (int(numParameters ** 0.5), int(numParameters ** 0.5))))),
    cmap=plt.get_cmap('seismic'),
    interpolation='none',
    extent=[-0.5, int(numParameters ** 0.5) - 0.5, -0.5, int(numParameters ** 0.5) - 0.5], vmin=0,
    vmax=1.2e-3)
axarr[0].invert_yaxis()
axarr[0].set_title("Mean of model blocks")
image1 = axarr[1].imshow(
    (1.0 / 450.0) * np.flipud(
        np.transpose(np.reshape(variances, (int(numParameters ** 0.5), int(numParameters ** 0.5))))),
    cmap=plt.get_cmap('cubehelix'),
    interpolation='none',
    extent=[-0.5, int(numParameters ** 0.5) - 0.5, int(numParameters ** 0.5) - 0.5, -0.5], vmin=0,
    vmax=3e-9)
axarr[1].invert_yaxis()
axarr[1].set_title("Variance of model blocks")

divider = make_axes_locatable(axarr[0])
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar0 = plt.colorbar(image0, cax=cax, format='%.1e')
cbar0.set_label('mean slowness[s/m]')

divider = make_axes_locatable(axarr[1])
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(image1, cax=cax, format='%.1e')
cbar1.set_label('variance')
# f.tight_layout()
# plt.tight_layout()
plt.tight_layout(pad=0, w_pad=0, h_pad=0)
# plt.subplots_adjust(wspace=0.5, hspace=-20)
if (sys.argv[2] == 0):
    plt.show()
else:
    print("Saving to", sys.argv[3])
    plt.savefig(sys.argv[3])
    plt.close()

