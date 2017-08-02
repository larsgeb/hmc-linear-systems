#!/usr/bin/env python
import numpy as np
# import random as random
import matplotlib.pyplot as plt

# ============================================================
# - Read data.
# ============================================================

fid = open('DATA/synthetics.txt')
data = fid.read().strip().split('\n')
fid.close()
fid = open('INPUT/setup.txt')
# setup = fid.read().strip().split('\n')
lines = fid.readlines()
setup = lines[1].split()
fid.close()

numRec = int(setup[0])
sepRec = float(setup[1])
firRec = float(setup[2])
travelTime = []
location = []

for x in range(0, data.__len__()):
    current = data[x].split(' ')
    travelTime.append(float(current[1]))
    location.append(float(current[0]))

# ============================================================
# - Plot data.
# ============================================================

# f, (ax1) = plt.subplots(1, 1, sharey=False)
ax1 = plt.gca()
plt.scatter(travelTime, location,marker=".")
ax1.set_title('travel time versus depth')
ax1.set_xlabel('time [s]')
ax1.set_ylabel('depth [m]')
plt.ylim([0-0.5*sepRec,location[location.__len__()-1]+0.5*sepRec])
plt.xlim([0,travelTime[travelTime.__len__()-1]])
ax1.invert_yaxis()
plt.show()
fid.close()