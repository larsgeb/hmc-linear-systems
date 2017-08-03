import numpy as np
import random as random
import matplotlib.pylab as pylab
import matplotlib.pyplot as plt

params = {'legend.fontsize': 'x-large',
          'figure.figsize': (6, 6),
          'axes.labelsize': 16,
          'axes.titlesize': 'x-large',
          'xtick.labelsize': 16,
          'ytick.labelsize': 16}
pylab.rcParams.update(params)

# ============================================================
# - Read samples and plot trajectory.
# ============================================================

fid = open('OUTPUT/multivariate.txt')
dummy = fid.read().strip().split()
fid.close()

x = []
y = []
for i in range(0, dummy.__len__() / 2):
    x.append(float(dummy[i * 2]))
    y.append(float(dummy[i * 2 + 1]))

plt.hist2d(x, y, bins=100, normed=True, cmap='binary')
plt.axis('equal')

plt.show()

mean_x = np.mean(x)
mean_y = np.mean(y)
cov_xx = 0.0
cov_yy = 0.0
cov_xy = 0.0

for i in range(x.__len__()):
    cov_xx += (mean_x - x[i]) * (mean_x - x[i])
    cov_yy += (mean_y - y[i]) * (mean_y - y[i])
    cov_xy += (mean_x - x[i]) * (mean_y - y[i])

cov_xx = cov_xx / (x.__len__())
cov_yy = cov_yy / (x.__len__())
cov_xy = cov_xy / (x.__len__())
print '---------------------------------------------------------------'
print 'mean_x=', mean_x, 'mean_y=', mean_y
print 'std_x=', np.sqrt(cov_xx), 'std_y=', np.sqrt(cov_yy), 'cov_xy=', cov_xy
print 'var_x=', (cov_xx), 'var_y=', (cov_yy)
print '---------------------------------------------------------------'
