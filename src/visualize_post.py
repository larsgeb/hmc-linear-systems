import matplotlib.pylab as pylab
import numpy as np
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

samples = np.loadtxt("OUTPUT/inversion_1.txt")

burnin = 1000

samples = samples[burnin:, :]

for i in range(0, 17):
    plt.plot(samples[:, i])

plt.show()

plt.plot(samples[:,-1])
plt.show()