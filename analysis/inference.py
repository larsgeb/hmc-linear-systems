# Visualization of Markov Chain Results
# Lars Gebraad, 2018

import sys
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

name = sys.argv[1]
dir_path = os.path.dirname(os.path.realpath(name))
samples = np.loadtxt(name)
burnin = 100
samples = samples[burnin::, :]
data = {}
means = []
names = []

zfill = int(np.ceil(np.log10(samples.shape[1] - 1)))

for i in range(0, samples.shape[1] - 1):
    means.append(np.mean(samples[:, i]))
    names.append("par" + str(i + 1).zfill(zfill))
    data["par" + str(i + 1).zfill(zfill)] = samples[:, i]

means = np.array(means)
means.shape = (means.shape[0], 1)

sns.set_style("darkgrid")
data = pd.DataFrame(data=data)

# Compute the correlation and covariance matrices
cor = data.corr()
cov = data.cov()
variances = np.diag(cov.values)
variances.shape = (variances.shape[0], 1)

# Save these matrices
np.savetxt(dir_path + "/correlation.txt", cor.values)
np.savetxt(dir_path + "/covariance.txt", cov.values)
np.savetxt(dir_path + "/means.txt", means)
np.savetxt(dir_path + "/variances.txt", variances)

# Custom cmap
cmap = sns.diverging_palette(220, 10, as_cmap=True)

# Plot both matrices
plt.subplots(figsize=(11, 9))
sns.heatmap(cor, cmap=cmap, vmax=1, vmin=-1, center=0, square=True, linewidths=.5, cbar_kws={"shrink": .5})
plt.show()
plt.subplots(figsize=(11, 9))
sns.heatmap(cov, cmap=cmap, vmax=1, vmin=-1, center=0, square=True, linewidths=.5, cbar_kws={"shrink": .5})
plt.show()

print("means:\t\t\t variances:")
print(np.hstack((means, variances)))

print("\r\n--------------------------------------------\r\n\r\n")
print("Correlation matrix:\r\n", cor)
print("\r\n--------------------------------------------\r\n\r\n")
print("Covariance matrix:\r\n", cov)

sns.set(style="ticks")
f, ax = plt.subplots(figsize=(7, 6))
plt.xscale('symlog', linthreshx=0.1)
plt.boxplot(samples[:, :-1], vert=False)
ax.xaxis.grid(True)
ax.set(ylabel="")
sns.despine(trim=True, left=True)
plt.xlabel('value (very skewed representation due to symlog axis)')
plt.yticks(np.arange(1, samples.shape[1]), names)
plt.show()
