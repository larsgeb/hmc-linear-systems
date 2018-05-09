import matplotlib.pylab as pylab
import numpy as np
import sys
import matplotlib.pyplot as plt

params = {'legend.fontsize': 'x-large',
          'figure.figsize': (5, 4),
          'axes.labelsize': 16,
          'axes.titlesize': 'x-large',
          'xtick.labelsize': 16,
          'ytick.labelsize': 16}
pylab.rcParams.update(params)

G = np.loadtxt("model_matrix.txt")

plt.imshow(G)
plt.colorbar()
plt.title("Forward model matrix")
plt.show()