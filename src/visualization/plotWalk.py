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



# If there's other stuff on the line, like misfit values
trailingElements = 1
nbi = 500
fid = open('OUTPUT/inversion_1.txt')
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

# Dimensions to visualize
dimA = 1
dimB = 2
with open('OUTPUT/inversion_1.txt') as f:
    iterations = sum(1 for _ in f) - nbi - 2
dimension = 4

x_plot = np.zeros(iterations)
y_plot = np.zeros(iterations)
x = np.zeros(iterations)
y = np.zeros(iterations)
qs = []

q_opt = np.zeros(dimension)
chi = 1.0e100

for i in range(0,iterations):

    x[i] = float(dummy[1 + dimA + (i + nbi) * (dimension + 1)])
    y[i] = float(dummy[1 + dimB + (i + nbi) * (dimension + 1)])
    x_plot[i] = x[i]
    y_plot[i] = y[i]

    # for parameter in range(0, dimension):
    #
    #     qs[parameter].append(float(dummy[2 + parameter + (i + nbi) * (dimension + 1)]))

    chi_test = float(dummy[2 + dimension + (i + nbi) * (dimension + 1)])
    if chi_test < chi:
        chi = chi_test
        print 'chi_min=', chi_test
        for k in range(dimension):
            q_opt[k] = float(dummy[2 + k + (i + nbi) * (dimension + 1)])

    if i > 0 and x[i] == x[i - 1] and x[i] > 0 and y[i] == y[i - 1]:
        x_plot[i] += epsilon_1 * random.gauss(0.0, 1.0)
        y_plot[i] += epsilon_2 * random.gauss(0.0, 1.0)

plt.plot(x_plot, y_plot, 'k', linewidth=0.05)
plt.plot(x_plot, y_plot, 'ro', linewidth=0.05, markersize=0.5)
axes = plt.gca()
plt.xlabel('parameter ' + str(dimA))
plt.ylabel('parameter ' + str(dimB))
plt.title('random walk')
plt.gcf().subplots_adjust(bottom=0.15)
plt.tight_layout()
# plt.savefig('OUTPUT/randomWalk.png')
# )
plt.show()
# plt.close()
# ============================================================
# - Histograms.
# ============================================================

xlimu = np.max(x)
xliml = np.min(x)
ylimu = np.max(y)
yliml = np.min(y)
plt.hist(x, bins=20, color='k', normed=True)
plt.xlim([xliml, xlimu])
plt.xlabel('m' + str(dimA))
plt.ylabel('posterior marginal')
plt.tight_layout()
plt.savefig('OUTPUT/marginal1.png')
plt.close()
# plt.show()
plt.hist(y, bins=20, color='k', normed=True)
plt.xlim([yliml, ylimu])
plt.xlabel('m' + str(dimB))
plt.ylabel('posterior marginal')
plt.tight_layout()
plt.savefig('OUTPUT/marginal2.png')
plt.close()
# plt.show()
plt.hist2d(x, y, bins=20, normed=True, cmap='binary')
# plt.axis('equal')
# plt.xlim([-20,40])
# plt.ylim([10,70])
plt.xlabel('m' + str(dimA))
plt.ylabel('m' + str(dimB))
plt.title('2D posterior marginal')
plt.colorbar()
plt.tight_layout()
plt.savefig('OUTPUT/marginal_2D.png')
plt.close()
# plt.show()
# ============================================================
# - Assess convergence.
# ============================================================
n = range(10, iterations, 10)

hist_final, bin_hist = np.histogram(x, bins=40, density=True)
diff = np.zeros(len(n))
k = 0
for i in n:
    hist, bin = np.histogram(x[0:i], bins=40, density=True)
    diff[k] = np.sqrt(np.sum((hist - hist_final) ** 2.0))
    k = k + 1

plt.semilogy(n, diff, 'k')
plt.xlabel('samples')
plt.ylabel('difference to final')
plt.tight_layout()
plt.savefig('OUTPUT/convergence1.png')
plt.close()

hist_final, bin = np.histogram(y, bins=40, density=True)
diff = np.zeros(len(n))

k = 0
for i in n:
    hist, bin = np.histogram(y[0:i], bins=40, density=True)
    diff[k] = np.sqrt(np.sum((hist - hist_final) ** 2.0))
    k = k + 1

plt.semilogy(n, diff, 'k')
plt.xlabel('samples')
plt.ylabel('difference to final')
plt.tight_layout()
plt.savefig('OUTPUT/convergence2.png')
plt.close()

## Last trajectory

