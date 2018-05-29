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

for i in range(0, samples.shape[1] - 1):
    data["par" + str(i + 1).zfill(3)] = samples[:, i]

sns.jointplot(x="par" + str(par_a).zfill(3), y="par" + str(par_b).zfill(3), data=data, kind="hex", color="#ffdaa8") \
    .set_axis_labels("par" + str(par_a).zfill(3), "par" + str(par_b).zfill(3))

plt.tight_layout()
plt.show()
