# Visualization of Markov Chain Joint Marginals
# Lars Gebraad, 2018

import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

name = sys.argv[1]
par_a = int(sys.argv[2])
par_b = int(sys.argv[3])

samples = np.loadtxt(name + "/samples.txt")
burnin = 100
samples = samples[burnin::, :]
data = {}
means = []

zfill = int(np.ceil(np.log10(samples.shape[1] - 1)))

for i in range(0, samples.shape[1] - 1):
    data["par" + str(i + 1).zfill(zfill)] = samples[:, i]

sns.jointplot(x="par" + str(par_a).zfill(zfill), y="par" + str(par_b).zfill(zfill), data=data, kind="hex", color="#ffdaa8") \
    .set_axis_labels("par" + str(par_a).zfill(zfill), "par" + str(par_b).zfill(zfill))

plt.tight_layout()
plt.show()
