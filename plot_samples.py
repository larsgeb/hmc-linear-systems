import numpy as np
import random as random
import matplotlib.pylab as pylab
import matplotlib.pyplot as plt

# plt.xkcd()
# ============================================================
# - Setup.
# ============================================================

# - Number of burn-in samples to be ignored.
nbi = 10
# - Dimensions of interest.
dim_1 = 0
dim_2 = 3
# - Incremental displacement for duplicate points.
epsilon_1 = 0.0003
epsilon_2 = 0.0003

params = {'legend.fontsize': 'x-large',
          'figure.figsize': (8, 8),
          'axes.labelsize': 20,
          'axes.titlesize': 'x-large',
          'xtick.labelsize': 20,
          'ytick.labelsize': 20}
pylab.rcParams.update(params)

# ============================================================
# - Read samples and plot trajectory.
# ============================================================

fid = open('OUTPUT/samples.txt')
dummy = fid.read().strip().split()
fid.close()
print
dimension = int(dummy[0])
iterations = int(dummy[dummy.__len__() - 1]) - nbi

x = np.zeros(iterations)
y = np.zeros(iterations)

x_plot = np.zeros(iterations)
y_plot = np.zeros(iterations)

q_opt = np.zeros(iterations)
chi = 1.0e100

for i in range(iterations):

    x[i] = float(dummy[2 + dim_1 + (i + nbi) * (dimension + 1)])
    y[i] = float(dummy[2 + dim_2 + (i + nbi) * (dimension + 1)])
    x_plot[i] = x[i];
    y_plot[i] = y[i];

    chi_test = float(dummy[2 + (dimension) + (i + nbi) * (dimension + 1)])
    if (chi_test < chi):
        chi = chi_test
        print 'chi_min=', chi_test
        for k in range(dimension):
            q_opt[k] = float(dummy[2 + k + (i + nbi) * (dimension + 1)])

    if (i > 0 and x[i] == x[i - 1] and x[i] > 0 and y[i] == y[i - 1]):
        x_plot[i] += epsilon_1 * random.gauss(0.0, 1.0)
        y_plot[i] += epsilon_2 * random.gauss(0.0, 1.0)

plt.plot(x_plot, y_plot, 'k', linewidth=0.05)
plt.plot(x_plot, y_plot, 'ro', linewidth=0.05, markersize=0.5)
# plt.axis('equal')
axes = plt.gca()
# axes.set_xlim([2315,2316])
# axes.set_ylim([495,496])
plt.xlabel('m' + str(dim_1 + 1))
plt.ylabel('m' + str(dim_2 + 1))
plt.title('random walk trajectory')
plt.savefig('OUTPUT/trajectory.png')
plt.savefig('OUTPUT/trajectory.pdf', format='pdf')
plt.close()
# plt.show()

# ============================================================
# - Histograms.
# ============================================================

xlimu = np.max(np.abs(x));
xliml = np.min(np.abs(x));
ylimu = np.max(np.abs(y));
yliml = np.min(np.abs(y));

plt.hist(x, bins=20, color='k', normed=True)
plt.xlim([xliml, xlimu])
plt.xlabel('m' + str(dim_1 + 1))
plt.ylabel('posterior marginal')
plt.savefig('OUTPUT/marginal1.png')
plt.close()
# plt.show()

plt.hist(y, bins=20, color='k', normed=True)
plt.xlim([yliml, ylimu])
plt.xlabel('m' + str(dim_2 + 1))
plt.ylabel('posterior marginal')
plt.savefig('OUTPUT/marginal2.png')
plt.close()
# plt.show()

plt.hist2d(x, y, bins=20, normed=True, cmap='binary')
# plt.axis('equal')
plt.xlabel('m' + str(dim_1 + 1))
plt.ylabel('m' + str(dim_2 + 1))
plt.title('2D posterior marginal')
plt.colorbar()
plt.savefig('OUTPUT/marginal_2D.png')
plt.close()
# plt.show()

# ============================================================
# - Assess convergence.
# ============================================================

n = range(10, iterations, 10)

hist_final, bin = np.histogram(x, bins=20, density=True)
diff = np.zeros(len(n))

k = 0
for i in n:
    hist, bin = np.histogram(x[0:i], bins=20, density=True)
    diff[k] = np.sqrt(np.sum((hist - hist_final) ** 2.0))
    k = k + 1

plt.semilogy(n, diff, 'k')
plt.xlabel('samples')
plt.ylabel('difference to final')
plt.savefig('OUTPUT/convergence1.png')
plt.close()

hist_final, bin = np.histogram(y, bins=20, density=True)
diff = np.zeros(len(n))

k = 0
for i in n:
    hist, bin = np.histogram(y[0:i], bins=20, density=True)
    diff[k] = np.sqrt(np.sum((hist - hist_final) ** 2.0))
    k = k + 1

plt.semilogy(n, diff, 'k')
plt.xlabel('samples')
plt.ylabel('difference to final')
plt.savefig('OUTPUT/convergence2.png')
plt.close()

# ============================================================
# - Make statistics.
# ============================================================

mean_x = np.mean(x)
mean_y = np.mean(y)

cov_xx = 0.0
cov_yy = 0.0
cov_xy = 0.0

for i in range(iterations - nbi):
    cov_xx += (mean_x - x[i]) * (mean_x - x[i])
    cov_yy += (mean_y - y[i]) * (mean_y - y[i])
    cov_xy += (mean_x - x[i]) * (mean_y - y[i])

cov_xx = cov_xx / (iterations)
cov_yy = cov_yy / (iterations)
cov_xy = cov_xy / (iterations)

print 'mean_x=', mean_x, 'mean_y=', mean_y
print 'std_xx=', np.sqrt(cov_xx), 'std_yy=', np.sqrt(cov_yy), 'cov_xy=', cov_xy
