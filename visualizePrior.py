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

mean = 1.0 / 1500.0
std = 1.0 / 1000.0
slowness = np.linspace(mean - 4 * std, mean + 10 * std, num=1000)
pdf = np.divide(
    np.exp(- np.divide(np.square(slowness - mean), (2.0 * std * std))),
    math.sqrt(2 * math.pi * std * std))

# Two subplots, the axes array is 1-d
f, axarr = plt.subplots(2)
axarr[0].plot(slowness, pdf)
axarr[0].set_xlim([mean - 4 * std, mean + 4 * std])

slowness = np.concatenate((np.linspace(mean - 0.9 * std, mean - 10 * std, num=10000),
                           np.linspace(mean + 10 * std, mean - 0.66 * std, num=10000)))
pdf = np.divide(
    np.exp(- np.divide(np.square(slowness - mean), (2.0 * std * std))),
    math.sqrt(2 * math.pi * std * std))

axarr[1].plot(1.0 / slowness, np.multiply(pdf, np.multiply(slowness, slowness)))
axarr[1].set_xlim([-4000, 4000])
plt.show()
