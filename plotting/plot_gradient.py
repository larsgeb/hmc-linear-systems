import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
import numpy as np
from numpy import ma

params = {'legend.fontsize': 'x-large',
          'figure.figsize': (8, 8),
          'axes.labelsize': 20,
          'axes.titlesize': 'x-large',
          'xtick.labelsize': 20,
          'ytick.labelsize': 20}
pylab.rcParams.update(params)

fid = open('OUTPUT/gradient.txt')
dummy = fid.read().strip().split("\n")
fid.close()

# parse data
q1 = []
q2 = []
dq1 = []
dq2 = []
for line in dummy:
    splitty = line.split()
    q1.append(float(splitty[0]))
    q2.append(float(splitty[1]))
    dq1.append(float(splitty[2]))
    dq2.append(float(splitty[3]))

plt.figure()
plt.title('Gradient of q1 and q2')
Q = plt.quiver(q1, q2, dq1, dq2)
plt.gca().set_aspect('equal', adjustable='box')
plt.show()