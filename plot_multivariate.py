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
trailingElements = 1
fid = open('OUTPUT/samples.txt')
dummy = fid.read().strip().split()
fid.close()
numParameters = int(dummy[0])
numSamples = int(dummy.__len__() - 2) / (numParameters + 1)
parameters = []
for i in range(0, numParameters):
    parameters.append([])

for i in range(0, numSamples):
    for index, parameter in enumerate(parameters):
        parameter.append(float(dummy[2 + index + (i * (numParameters + trailingElements))]))

means = []
for index, parameter in enumerate(parameters):
    means.append(np.mean(parameter))

covariances = [[0] * numParameters for i in range(numParameters)]
for index1, parameter1 in enumerate(parameters):
    for index2, parameter2 in enumerate(parameters):
        covariances[index1][index2] = 0
        for indexSample in range(numSamples):
            covariances[index1][index2] += (means[index1] - parameter1[indexSample]) \
                                           * (means[index2] - parameter2[indexSample])
        covariances[index1][index2] = (covariances[index1][index2] / numSamples)

print 'Means;'
printMatrixE(np.matrix(means))
print 'Covariances;'
printMatrixE(np.matrix(covariances))
print "Standard deviations;"
for index, covariance in enumerate(covariances):
    print "Parameter", index + 1, np.sqrt(covariance[index])

plt.imshow(covariances, cmap=plt.get_cmap('binary'), interpolation='none',
           extent=[0.5, numParameters + 0.5, numParameters + 0.5, 0.5])
plt.show()
