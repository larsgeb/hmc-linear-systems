import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def autocorr(x):
    result = np.correlate(x, x, mode='full')
    return result[int(result.size / 2):] / (np.max(result))


samples = np.loadtxt("OUTPUT/inversion_1.txt")

burnin = 0

sns.set_style("darkgrid")
samples = samples[burnin::1000, :]

print(samples.shape)

data = {}

for i in range(0, 17):
    plt.plot(samples[:, i])
    data[str(i).zfill(2)] = samples[:, i]

plt.show()

plt.plot(samples[:, -1], 'black')
plt.xlabel('sample')
plt.ylabel('misfit')
plt.show()

data = pd.DataFrame(data=data)

# Compute the correlation matrix
corr = data.corr()

# Generate a mask for the upper triangle
mask = np.zeros_like(corr, dtype=np.bool)
mask[np.triu_indices_from(mask)] = True

# Set up the matplotlib figure
f, ax = plt.subplots(figsize=(11, 9))

# Generate a custom diverging colormap
cmap = sns.diverging_palette(220, 10, as_cmap=True)

# Draw the heatmap with the mask and correct aspect ratio
sns.heatmap(corr, mask=mask, cmap=cmap, vmax=1, vmin=-1,
            center=0, square=True, linewidths=.5, cbar_kws={"shrink": .5})

plt.show()

print(corr)

sns.set_style("darkgrid")
plt.subplot2grid((3,3), (0,0), colspan=2, rowspan=3)
for i in range(0, 17):
    au = autocorr(samples[:, i] - np.mean(samples[:, i]))
    plt.plot(au, color='black')
plt.xlabel("Sample lag (k)")
plt.ylabel("Corellation")
plt.subplot2grid((3,3), (0,2), colspan=1, rowspan=3)
for i in range(0, 17):
    au = autocorr(samples[:, i] - np.mean(samples[:, i]))
    plt.plot(au[0:40], color='black')
plt.xlabel("Sample lag (k)")
plt.show()
