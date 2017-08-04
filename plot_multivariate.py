import matplotlib.pylab as pylab
import numpy as np
import matplotlib.pyplot as plt


def printMatrixE(a):
    print "Matrix[" + ("%d" % a.shape[0]) + "][" + ("%d" % a.shape[1]) + "]"
    rows = a.shape[0]
    cols = a.shape[1]
    for i in range(0, rows):
        for j in range(0, cols):
            print("%6.10f" % a[i, j]),
        print
    print


params = {'legend.fontsize': 'x-large',
          'figure.figsize': (6, 6),
          'axes.labelsize': 16,
          'axes.titlesize': 'x-large',
          'xtick.labelsize': 16,
          'ytick.labelsize': 16}
pylab.rcParams.update(params)

# If there's other stuff on the line, like misfit values
trailingElements = 0
nbi = 0
fid = open('OUTPUT/multivariate.txt')
dummy = fid.read().strip().split()
fid.close()
numParameters = int(dummy[0])
numSamples = int(dummy.__len__() - 2) / (numParameters + 1)
parameters = []
for i in range(0, numParameters):
    parameters.append([])

for i in range(nbi, numSamples):
    for index, parameter in enumerate(parameters):
        parameter.append(float(dummy[2 + index + (i * (numParameters + trailingElements))]))

means = []
for index, parameter in enumerate(parameters):
    means.append(np.mean(parameter))

covariances = [[0] * numParameters for i in range(numParameters)]
for index1, parameter1 in enumerate(parameters):
    for index2, parameter2 in enumerate(parameters):
        covariances[index1][index2] = 0
        for indexSample in range(numSamples - nbi):
            covariances[index1][index2] += (means[index1] - parameter1[indexSample]) *\
                                           (means[index2] - parameter2[indexSample])
        covariances[index1][index2] = (covariances[index1][index2] / (numSamples - nbi))

# covariances_alt = [[0] * numParameters for i in range(numParameters)]
# for index1, parameter1 in enumerate(parameters):
#     for index2, parameter2 in enumerate(parameters):
#         covariances_alt[index1][index2] = 0
#         for indexSample1 in range(numSamples - nbi):
#             for indexSample2 in range(numSamples - nbi):
#                 covariances_alt[index1][index2] += (parameter1[indexSample2] - parameter1[indexSample1]) * (
#                     parameter2[indexSample2] - parameter2[indexSample1])
#         covariances_alt[index1][index2] = (covariances[index1][index2] / (numSamples - nbi) * (numSamples - nbi))

print 'Means;'
printMatrixE(np.matrix(means))
print 'Inverse means;'
printMatrixE(1 / np.matrix(means))
print 'Covariances;'
printMatrixE(np.matrix(covariances))
# print 'ALternative covariances;'
# printMatrixE(np.matrix(covariances_alt))
print "Standard deviations;"
for index, covariance in enumerate(covariances):
    print "Parameter", index + 1, np.sqrt(covariance[index])

maxCov = np.max(np.abs(covariances))

plt.imshow(covariances, cmap=plt.get_cmap('seismic'), interpolation='none',
           extent=[0.5, numParameters + 0.5, numParameters + 0.5, 0.5], vmin=-maxCov, vmax=maxCov)
plt.show()
