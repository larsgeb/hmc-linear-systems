import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

name = "inversion2"

samples = np.loadtxt(name + "/samples.txt")

burnin = 10

samples = samples[burnin::, :]

data = {}

means = []

for i in range(0, 17):
    means.append(np.mean(samples[:, i]))
    data["par" + str(i + 1).zfill(2)] = samples[:, i]

means = np.array(means)
means.shape=(means.shape[0],1)

sns.set_style("darkgrid")
print("Means:\r\n", (means), "\r\n")
data = pd.DataFrame(data=data)

# Compute the correlation matrix
corr = data.corr()

# Set up the matplotlib figure
f, ax = plt.subplots(figsize=(11, 9))

# Generate a custom diverging colormap
cmap = sns.diverging_palette(220, 10, as_cmap=True)

sns.heatmap(corr, cmap=cmap, vmax=1, vmin=-1,
            center=0, square=True, linewidths=.5, cbar_kws={"shrink": .5})

plt.show()

cov = data.cov()
# Set up the matplotlib figure
f, ax = plt.subplots(figsize=(11, 9))
# Generate a custom diverging colormap
cmap = sns.diverging_palette(220, 10, as_cmap=True)
sns.heatmap(cov, cmap=cmap, vmax=1, vmin=-1,
            center=0, square=True, linewidths=.5, cbar_kws={"shrink": .5})
plt.show()

print("\r\n--------------------------------------------\r\n\r\n")
print("Correlation matrix:\r\n",corr)
print("\r\n--------------------------------------------\r\n\r\n")
print("Covariance matrix:\r\n", cov)
