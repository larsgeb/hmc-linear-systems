import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

name = sys.argv[1]
samples = np.loadtxt(name + "/samples.txt")
burnin = 100
samples = samples[burnin::, :]
data = {}
means = []

for i in range(0, 17):
    means.append(np.mean(samples[:, i]))
    data["par" + str(i + 1).zfill(2)] = samples[:, i]

means = np.array(means)
means.shape = (means.shape[0], 1)

sns.set_style("darkgrid")
data = pd.DataFrame(data=data)

# Compute the correlation and covariance matrices
cor = data.corr()
cov = data.cov()
variances = np.diag(cov.as_matrix())
variances.shape = (variances.shape[0], 1)

# Save these matrices
np.savetxt(name + "/correlation.txt", cor.as_matrix())
np.savetxt(name + "/covariance.txt", cov.as_matrix())

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
