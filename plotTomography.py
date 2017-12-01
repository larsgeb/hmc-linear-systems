import matplotlib.pylab as pylab
import numpy as np
import matplotlib.pyplot as plt
import math
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec

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
fid = open('OUTPUT/samples_neal_other_tuning.txt')
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
# plt.imshow(np.transpose(np.reshape(variances, (101,101))), cmap=plt.get_cmap('seismic'), interpolation='none',
#            extent=[-0.5, 10 + 0.5, 10 + 0.5, -0.5], vmin=0, vmax=maxVar)
# cbar = plt.colorbar(ticks=np.arange(0, maxVar, maxVar/10))
# # plt.gca().invert_yaxis()
# plt.savefig("OUTPUT/tomographyVariances.pdf", format='pdf')
# plt.close()
#
# plt.imshow(1 / np.transpose(np.reshape(means, (101,101))), cmap=plt.get_cmap('binary'), interpolation='none',
#            extent=[-0.5, 10 + 0.5, -0.5, 10 + 0.5], vmin=0000, vmax=3001)
# cbar = plt.colorbar(ticks=np.arange(0, 3001, 500))
# plt.gca().invert_yaxis()
# cbar.set_label('Speed of sound [m/s]')
# plt.savefig("OUTPUT/tomographyMeans.pdf", format='pdf')
# plt.close()

f, axarr = plt.subplots(1, 2)
# gs1 = gridspec.GridSpec(1, 2)
# gs1.update(wspace=50, hspace=50)
image0 = axarr[0].imshow(np.transpose(np.reshape(means, (11, 11))), cmap=plt.get_cmap('binary'),
                         interpolation='none',
                         extent=[-0.5, 101 + 0.5, -0.5, 101 + 0.5])
# image0 = axarr[0].imshow(1 / np.transpose(np.reshape(means, (101, 101))), cmap=plt.get_cmap('binary'),
#                          interpolation='none',
#                          extent=[-0.5, 101 + 0.5, -0.5, 101 + 0.5], vmin=0000, vmax=3001)
axarr[0].invert_yaxis()
axarr[0].set_title("Mean of model blocks")
image1 = axarr[1].imshow(np.transpose(np.reshape(variances, (11, 11))), cmap=plt.get_cmap('seismic'),
                         interpolation='none',
                         extent=[-0.5, 101 + 0.5, 101 + 0.5, -0.5], vmin=0, vmax=maxVar)
axarr[1].invert_yaxis()
axarr[1].set_title("Variance of model blocks")
f.suptitle("Inversion of tomographic data. Underlying model: checkerboard of 1000/2000 m/s cells\n10201 parameters, "
           "100 samples (all accepted), temperature=10,000")


divider = make_axes_locatable(axarr[0])
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar0 = plt.colorbar(image0, cax=cax)
cbar0.set_label('Slowness of sound [s/m]')

divider = make_axes_locatable(axarr[1])
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(image1, cax=cax)
cbar1.set_label('Variance of q')
# f.tight_layout()

plt.subplots_adjust(wspace=0.5, hspace=20)
plt.show()
