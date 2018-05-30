# Visualization of Markov Chain Behaviour
# Lars Gebraad, 2018

import sys
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Find next power of two from any given number
def next_pow_two(x):
    return 1<<(x-1).bit_length()



# Rather fast implementation of autocorrelation using FFT
def autocorrFFT(x, norm=True):
    x = np.atleast_1d(x)
    if len(x.shape) != 1:
        raise ValueError("invalid dimensions for 1D autocorrelation function")
    n = next_pow_two(len(x))

    # Compute the FFT and then (from that) the auto-correlation function
    f = np.fft.fft(x - np.mean(x), n=2 * n)
    acf = np.fft.ifft(f * np.conjugate(f))[:len(x)].real
    acf /= 4 * n

    # Optionally normalize
    if norm:
        acf /= acf[0]

    return acf

name = sys.argv[1]
samples = np.loadtxt(name + "/samples.txt")
burnin = 100
samples = samples[burnin::, :]

zfill = int(np.ceil(np.log10(samples.shape[1] - 1)))

# Samples plot
legend = []
sns.set_style("darkgrid")
for i in range(0, 17):
    plt.plot(samples[:, i])
    legend.append("par" + str(i + 1).zfill(zfill))
legend = plt.legend(legend)
plt.show()

# Sample misfit plot
sns.set_style("darkgrid")
plt.plot(samples[:, -1], 'black')
plt.xlabel('sample')
plt.ylabel('misfit')
plt.show()

# Sample autocorrelation plot
sns.set_style("darkgrid")
plt.subplot2grid((3, 3), (0, 0), colspan=2, rowspan=3)
for i in range(0, 17):
    au = autocorrFFT(samples[::, i] - np.mean(samples[:, i]))
    plt.plot(au, color='black')
plt.xlabel("Sample lag (k)")
plt.ylabel("autocorrelation")
plt.subplot2grid((3, 3), (0, 2), colspan=1, rowspan=3)
for i in range(0, 17):  # One could save some time here by saving the result of the previous loop
    au = autocorrFFT(samples[::, i] - np.mean(samples[:, i]))
    plt.plot(au[0:40], color='black')
plt.xlabel("Sample lag (k)")
plt.show()
