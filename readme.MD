

    Services
    Documents
        Untitled Document.md

    Preview as
    Export as
    Save to
    Import from

Document Name

Words: 246
Markdown
Toggle Zen Mode
Preview
HMC Inversion of simple linear models

This set of programs is built to incorporate MCMC inversion using Hamiltonian Mechanics. Essentially, any amount of
parameters and any linear forward model matrix can be supplied, as long as one can describe the measurement errors, data
uncertainty and misfit function. This could be used as a framework for building any inversion using HMC.
Compilation

All build targets can easily be prepared with CMake. Using at least CMake 3.7 allows you to run (from the
root directory);

cmake .

This prepares the platform dependent compilation instructions. For example, on Linux one would now run for
all separate targets:

make <target>

There are several build targets:

    generateSynthetics.
    This program will read the forward matrix in INPUT/matrix.txt and generate synthetics which are saved in synthetics.txt. Model parameters are as of yet hardcoded.

    runSampling.
    This program will take samples from the posterior distribution described by hardcoded prior information, the forward model from INPUT/matrix.txt and synthetic data INPUT/synthetics.txt (with hardcoded data-covariance). The output are the samples in OUTPUT/samples.txt and the trajectory of the last sample in OUTPUT/trajectory.txt. Some diagnostics are written to the console.

    generateGradient.
    This programs evaluates misfit gradients on a hardcoded interval, using the same input as runSampling.

    testMultivariate.
    This program draws from hardcoded multivariate correlated gaussians, writing it’s samples to OUTPUT/multivariate.txt. Used to evaluate the validity of the Cholesky decomposition.

Code documentation

Documentation generated using doxygen is found in th docs/ subdirectory, but is also found hosted at my GitHub Pages.
